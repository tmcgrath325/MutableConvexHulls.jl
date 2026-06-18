# jarvissortedsearch is a potential optimization of jarvissearch that exploits the
# sorted order of the hull list. It is not currently used (jarvissearch is used instead)
# but is retained as a reference for future optimization work.
#
# Precondition: pointslist must be a strictly-convex partial (lower or upper) hull list,
# i.e. built with collinear=false. For such lists, the betterturn quality function is
# unimodal over the sorted vertices: it rises to a single peak then falls. This function
# exploits that unimodality to exit early once both the predecessor and successor are
# worse than the current node, avoiding a full scan of the list.
#
# The unimodality assumption fails for collinear=true hull lists (many collinear points
# make the quality non-monotone) and for full convex hull lists (the circular structure
# and non-transitivity of betterturn can produce cyclic comparisons).
function jarvissortedsearch(query::AbstractNode, prevedge, pointslist::AbstractLinkedList, betterturn::Function)
    length(pointslist) == 0 && throw(ArgumentError("The list of points must not be empty."))
    length(pointslist) == 1 && return head(pointslist)
    # Find the first and last nodes with coordinates distinct from the query.
    # prevnode (from the tail) provides the circular-boundary comparison for the head.
    prevnode = getfirst(x -> !coordsareequal(query.data, x.data), ListNodeIterator(pointslist; rev=true))
    prevnode === nothing && return head(pointslist)
    currentnode = getfirst(x -> !coordsareequal(query.data, x.data), ListNodeIterator(pointslist; rev=false))
    currentnode === prevnode && return currentnode
    # Comparing currentnode (head) against prevnode (tail) handles the wrap-around case:
    # if the peak of the unimodal function sits at the very start of the list, both the
    # tail (prev) and the second node (next) will be worse, and we return immediately.
    prev_worse = !betterturn(prevedge, query.data, currentnode.data, prevnode.data)
    nextnode = getfirst(x -> !coordsareequal(query.data, x.data), ListNodeIterator(currentnode.next; rev=false))
    while nextnode !== nothing
        next_worse = !betterturn(prevedge, query.data, currentnode.data, nextnode.data)
        if prev_worse && next_worse
            return currentnode
        end
        prev_worse = !next_worse
        prevnode = currentnode
        currentnode = nextnode
        nextnode = getfirst(x -> !coordsareequal(query.data, x.data), ListNodeIterator(currentnode.next; rev=false))
    end
    return currentnode
end

# Dispatch to jarvissortedsearch when the hull list is a strictly-convex partial (lower
# or upper) hull built with collinear=false; fall back to jarvissearch otherwise.
function hulllistsearch(query, prevedge, htarget::AbstractLinkedList, betterturn, partial::Bool, collinear::Bool)
    partial && !collinear ?
        jarvissortedsearch(query, prevedge, htarget, betterturn) :
        jarvissearch(query, prevedge, ListNodeIterator(htarget), betterturn)
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

    targetscollinear = [hl.collinear for hl in hulls]

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
    merge_hull_lists!(mergedhull, hulltargets, buildinreverse(h), h.orientation, h.collinear, h.sortedby, targetscollinear, partial, upper)
    return h
end
"""
    mergehulls(hull, otherhulls...)

Return a new hull of the same type as `hull` containing the points of `hull` and
`otherhulls`, without mutating any argument. See [`mergehulls!`](@ref) for the
in-place form.
"""
mergehulls(h::H, others::H...) where {H <: Union{MutableConvexHull, MutableLowerConvexHull, MutableUpperConvexHull}} = mergehulls!(copy(h), others...)

function merge_hull_lists!(mergedhull::AbstractList, hulltargets::Vector{<:AbstractList}, rev::Bool, orientation::HullOrientation, collinear::Bool, sortedby::Function, targetscollinear::Vector{Bool}, partial::Bool, upper::Bool)
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

    # perform jarvis march with search that makes use of the sorted nature of the hulls
    counter = 0
    current = start
    currentdata = current.data
    prevdata = currentdata
    candidates = [head(x) for x in hulltargets]
    prevedge = upper ? UP : DOWN
    while counter == 0 || current !== stop
        if counter > maxlength
            throw(ErrorException("More points were added to the hull ($counter) than exist in the original hulls to be merged ($maxlength)."))
        end
        counter += 1
        for (i, ht) in enumerate(hulltargets)
            candidates[i] = hulllistsearch(current, prevedge, ht, betterturn, partial, targetscollinear[i])
        end
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

    # if these are only partial convex hulls (i.e. upper or lower), add the stopping point at the end of the hull
    if partial
        push!(mergedhull, stop.data)
        addtarget!(tail(mergedhull), stop.target)
    end
    return mergedhull
end

function merge_hull_lists!(h::AbstractChanConvexHull)
    mergedhull = h.hull
    hulltargets = filter(!isempty, [hl.hull for hl in h.subhulls])
    rev = buildinreverse(h.subhulls[1])
    targetscollinear = fill(h.collinear, length(h.subhulls))
    upper = eltype(h.subhulls) <: MutableUpperConvexHull
    partial = upper || eltype(h.subhulls) <: MutableLowerConvexHull
    merge_hull_lists!(mergedhull, hulltargets, rev, h.orientation, h.collinear, h.sortedby, targetscollinear, partial, upper)
    return h
end

function fallback_merge_hull_lists!(mergedhull::AbstractList, hulltargets::Vector{<:AbstractList}, rev::Bool, orientation::HullOrientation, collinear::Bool, sortedby::Function, targetscollinear::Vector{Bool}, partial::Bool, upper::Bool)
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
    targetscollinear = fill(h.collinear, length(h.subhulls))
    upper = eltype(h.subhulls) <: MutableUpperConvexHull
    partial = upper || eltype(h.subhulls) <: MutableLowerConvexHull
    fallback_merge_hull_lists!(mergedhull, hulltargets, rev, h.orientation, h.collinear, h.sortedby, targetscollinear, partial, upper)
    return h
end
