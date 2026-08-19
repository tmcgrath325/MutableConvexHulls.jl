# The first node of `htarget` lying ahead of `key`, or its tail if none does.
#
# A partial hull's vertices ascend in `sortedby` (descend when built in
# reverse), so only the nodes beyond the march's current key can be the next
# one. Restricting the search to them is what keeps the march from doubling
# back: where two candidates lie closer together than the arithmetic can
# resolve — coordinates spanning many orders of magnitude leave the cross
# product below the precision of even the extended-precision fallback — the
# turn predicate reads them as collinear, and its tie-break among collinear
# candidates prefers the more distant of the two, which may sit behind the
# march.
#
# The march's key only ever advances, so each list is scanned from where its
# last search left off and every node is passed over at most once across the
# whole march.
function aheadfrom(cursor::AbstractNode, sortedby::Function, key, rev::Bool)
    while !attail(cursor) && !(rev ? sortedby(cursor.data) < key : sortedby(cursor.data) > key)
        cursor = cursor.next
    end
    return cursor
end

"""
    mergehulls!(hull, otherhulls...)

Merge the points contained in `otherhulls` into `hull` and return `hull`. All
arguments must be the same concrete hull type (`MutableConvexHull`,
`MutableLowerConvexHull`, or `MutableUpperConvexHull`); [`AbstractChanConvexHull`](@ref)
subtypes are not accepted. See [Chan's algorithm](https://en.wikipedia.org/wiki/Chan%27s_algorithm)
for a similar approach.
"""
function mergehulls!(h::H, others::H...) where {H <: Union{MutableConvexHull, MutableLowerConvexHull, MutableUpperConvexHull}}
    mergedhull = h.hull
    mergedpoints = h.points

    # filter out empty hulls
    hulls = filter(x -> length(x.hull) > 0, [h, others...])
    length(hulls) == 0 && return h

    # Set up copies of the hulls that point to the new points list
    hulltargets = [TargetedLinkedList(mergedpoints) for i in 1:length(hulls)]

    for (originalhull, htarget) in zip(hulls, hulltargets)
        for hullnode in ListNodeIterator(originalhull.hull)
            push!(htarget, hullnode.data)
            tail(htarget).target = hullnode.target
        end
    end

    # add points from all hulls into the points list.
    for originalhull in hulls
        if originalhull !== h
            for pointnode in ListNodeIterator(originalhull.points)
                removetarget!(pointnode)
                pointnode.list = mergedpoints
                pointnode.up = pointnode
                pointnode.down = pointnode
                push!(mergedpoints, pointnode)
            end
        end
    end

    # merge the hull points into a new convex hull
    upper = H <: MutableUpperConvexHull
    partial = upper || H <: MutableLowerConvexHull
    merge_hull_lists!(mergedhull, hulltargets, buildinreverse(h), h.orientation, h.collinear, h.sortedby, partial, upper)
    return h
end
"""
    mergehulls(hull, otherhulls...)

Return a new hull of the same type as `hull` containing the points of `hull` and
`otherhulls`, without mutating any argument. See [`mergehulls!`](@ref) for the
in-place form.
"""
mergehulls(h::H, others::H...) where {H <: Union{MutableConvexHull, MutableLowerConvexHull, MutableUpperConvexHull}} = mergehulls!(copy(h), others...)

function merge_hull_lists!(mergedhull::AbstractList, hulltargets::Vector{<:AbstractList}, rev::Bool, orientation::HullOrientation, collinear::Bool, sortedby::Function, partial::Bool, upper::Bool)
    empty!(mergedhull) # start with an empty hull
    # handle simple cases
    isempty(hulltargets) && return empty!(mergedhull)
    if length(hulltargets) == 1
        ht = only(hulltargets)
        if (length(ht) == 1)
            hthead = head(ht)
            push!(mergedhull, hthead.data)
            addtarget!(tail(mergedhull), hthead.target)
            return mergedhull
        end
    end
    # determine starting and stopping points for general case
    f = x -> sortedby(x.data)
    start = rev ? argmax(f, [head(x) for x in hulltargets]) :
        argmin(f, [head(x) for x in hulltargets])
    stop = start
    if partial
        stop = rev ? argmin(f, [tail(x) for x in hulltargets]) :
            argmax(f, [tail(x) for x in hulltargets])
    end
    stopdata = stop.data

    maxlength = sum(length, hulltargets)

    # add first point to hull
    pushfirst!(mergedhull, start.data)
    addtarget!(head(mergedhull), start.target)

    maxlength <= 1 && return mergedhull

    # prepare orientation test
    betterturn(prevedge, o, a, b) = collinear ? iscloserturn(!orientation, prevedge, o, a, b) : isfurtherturn(!orientation, prevedge, o, a, b)

    # perform jarvis march
    counter = 0
    current = start
    currentdata = current.data
    prevdata = currentdata
    candidates = [head(x) for x in hulltargets]
    cursors = [head(x) for x in hulltargets]
    prevedge = upper ? UP : DOWN
    while counter == 0 || current !== stop
        if counter > maxlength
            throw(ErrorException("More points were added to the hull ($counter) than exist in the original hulls to be merged ($maxlength)."))
        end
        counter += 1
        empty!(candidates)
        key = sortedby(current.data)
        for (i, ht) in enumerate(hulltargets)
            if partial
                cursors[i] = aheadfrom(cursors[i], sortedby, key, rev)
                attail(cursors[i]) && continue          # this list holds nothing ahead
                push!(candidates, jarvissearch(current, prevedge, ListNodeIterator(cursors[i]), betterturn))
            else
                # A full hull wraps around, so its march revisits keys it has
                # passed and every vertex stays a candidate.
                push!(candidates, jarvissearch(current, prevedge, ListNodeIterator(ht), betterturn))
            end
        end
        isempty(candidates) && break        # nothing lies ahead: the chain is complete
        next = jarvissearch(current, prevedge, candidates, betterturn)
        if coordsareequal(current.data, next.data)
            if length(mergedhull) == 1
                return mergedhull
            else
                throw(ErrorException("Jarvis March failed to progress."))
            end
        end
        nextdata = next.data
        # stop adding points when the stopping point has been reached
        if coordsareequal(prevdata, nextdata) || coordsareequal(stopdata, nextdata)
            break
        end
        prevdata = currentdata
        currentdata = nextdata
        # add the next node to the hull
        push!(mergedhull, next.data)
        addtarget!(tail(mergedhull), next.target)
        prevedge = sub2d(next.data, current.data)
        current = next
    end

    # A partial hull (upper or lower) ends at its stopping point. The march
    # reaches it only when something lies beyond the vertex it is on, so where
    # every remaining point shares those coordinates the chain already ends
    # there and appending it again would put two vertices on one point.
    if partial && !coordsareequal(tail(mergedhull).data, stop.data)
        push!(mergedhull, stop.data)
        addtarget!(tail(mergedhull), stop.target)
    end
    return mergedhull
end

function merge_hull_lists!(h::AbstractChanConvexHull)
    mergedhull = h.hull
    hulltargets = filter(!isempty, [hl.hull for hl in h.subhulls])
    rev = buildinreverse(h.subhulls[1])
    upper = eltype(h.subhulls) <: MutableUpperConvexHull
    partial = upper || eltype(h.subhulls) <: MutableLowerConvexHull
    merge_hull_lists!(mergedhull, hulltargets, rev, h.orientation, h.collinear, h.sortedby, partial, upper)
    return h
end

function fallback_merge_hull_lists!(mergedhull::AbstractList, hulltargets::Vector{<:AbstractList}, rev::Bool, orientation::HullOrientation, collinear::Bool, sortedby::Function, partial::Bool, upper::Bool)
    empty!(mergedhull)
    all_points_nodes = sort(collect(Iterators.flatten([n.target for n in ListNodeIterator(h; rev = rev)] for h in hulltargets)); by = x -> sortedby(x.data), rev = rev)
    return if !isempty(all_points_nodes)
        stop = first(all_points_nodes)
        if partial
            stop = rev ? argmin(x -> sortedby(x.data), all_points_nodes) : argmax(x -> sortedby(x.data), all_points_nodes)
        end
        push!(mergedhull, first(all_points_nodes).data)
        addtarget!(head(mergedhull), first(all_points_nodes))
        jarvismarch!(mergedhull, all_points_nodes, collinear, orientation, upper ? UP : DOWN, stop)
        if partial && head(mergedhull).target !== stop
            push!(mergedhull, stop.data)
            addtarget!(tail(mergedhull), stop)
        end
    end
end

function fallback_merge_hull_lists!(h::AbstractChanConvexHull)
    mergedhull = h.hull
    hulltargets = filter(!isempty, [hl.hull for hl in h.subhulls])
    rev = buildinreverse(h.subhulls[1])
    upper = eltype(h.subhulls) <: MutableUpperConvexHull
    partial = upper || eltype(h.subhulls) <: MutableLowerConvexHull
    fallback_merge_hull_lists!(mergedhull, hulltargets, rev, h.orientation, h.collinear, h.sortedby, partial, upper)
    return h
end
