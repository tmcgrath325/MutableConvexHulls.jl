function jarvissortedsearch(query::AbstractNode, prevedge, pointslist::AbstractLinkedList, betterturn::Function)
    pointslist.len == 0 && throw(ArgumentError("The list of points must not be empty."))
    pointslist.len == 1 && return head(pointslist)
    head(pointslist).data == tail(pointslist).data && throw(ArgumentError("All points in the list are duplicates."))

    prev_worse = betterturn(prevedge, query.data, tail(pointslist).data, head(pointslist).data)
    for target in ListNodeIterator(pointslist)
        if query.data == target.next.data
            continue
        end
        next_worse = !betterturn(prevedge, query.data, target.data, target.next.data)
        if prev_worse && next_worse
            return target
        end
        prev_worse = !next_worse
    end
    return tail(pointslist)
end

# Note: binary search is a bad idea for linked lists

# function jarvisbinarysearch(query::PointNode, prevedge, pointslist::AbstractLinkedList, betterturn::Function)
#     pointslist.len == 0 && throw(ArgumentError("The list of points must not be empty."))
#     pointslist.len == 1 && return head(pointslist)
#     head(pointslist).data == tail(pointslist).data && throw(ArgumentError("All points in the list are duplicates."))

#     # initialize interval bounds
#     left = 1
#     right = pointslist.len
#     middle = Int(floor(pointslist.len/2))
#     target = getnode(pointslist, middle)

#     # recursively shrink interval until a tangent is found
#     return jarvisbinarysearchinterval(left, right, middle, target, query, prevedge, pointslist, betterturn)
# end

# # Perform binary search in the interval specified by left and right. The middle index and target node are already determined for the interval
# function jarvisbinarysearchinterval(left::Int, right::Int, middle::Int, target::AbstractNode, query::PointNode, prevedge, pointslist::AbstractLinkedList, betterturn::Function)
#     left == middle == right && return target

#     # wrap around the list to get prev/next values for the ends
#     prev = middle == 1 ? tail(pointslist) : target.prev
#     next = middle == pointslist.len ? head(pointslist) : target.next

#     # if the query is a duplicate of a hull point, 
#     (query.data == target.data) && return next
#     (query.data == prev.data) && return target
#     (query.data == next.data) && return next.next

#     # check if the adjacent nodes represent "better" options than the target one
#     prev_better = betterturn(prevedge, query.data, target.data, prev.data)
#     next_better = betterturn(prevedge, query.data, target.data, next.data)
#     if !prev_better && !next_better # the target is the best choice
#         return target
#     elseif right - left == 1
#         # handle cases when the best point lies "before" the beginning of the list
#         if middle == 1
#             while prev_better
#                 target = prev
#                 prev = target.prev
#                 prev_better = betterturn(prevedge, query.data, target.data, prev.data)
#             end
#             return target
#         end
#         # handle cases when the best point lies "after" the end of the list
#         if middle == pointslist.len
#             while next_better
#                 target = next
#                 next = target.next
#                 next_better = betterturn(prevedge, query.data, target.data, next.data)
#             end
#             return target
#         end
#         # otherwise, pick the better of the two remaining points
#         return prev_better ? prev : next
#     elseif prev_better # explore the left interval
#         (right - left) == 0 && return prev
#         right = middle
#         middle = Int(floor((right - left)/2) + left)
#         idx_change = right - middle
#         for i=1:idx_change
#             target = target.prev
#         end
#     else # explore the right interval
#         (right - left) == 0 && return next
#         left = middle
#         middle = Int(ceil((right - left)/2) + left)
#         idx_change = middle - left
#         for i=1:idx_change
#             target = target.next
#         end
#     end

#     return jarvisbinarysearchinterval(left, right, middle, target, query, prevedge, pointslist, betterturn)
# end

"""
    mergehulls!(hull, otherhulls...)

Merge the points contained in `otherhulls` into `hull`. See [Chan's algorithm](https://en.wikipedia.org/wiki/Chan%27s_algorithm)
for a similar approach.
"""
function mergehulls!(h::H, others::H...) where H<:AbstractConvexHull
    mergedhull = h.hull
    mergedpoints = h.points
    
    # filter out empty hulls
    hulls = filter(x->length(x.hull)>0,[h, others...])
    length(hulls) == 0 && return h

    targetscollinear = [hl.collinear for hl in hulls]

    # Set up copies of the hulls that point to the new points list
    hulltargets = [TargetedLinkedList(mergedpoints) for i=1:length(hulls)]
    for (originalhull,htarget) in zip(hulls,hulltargets)
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
mergehulls(h::H, others::H...) where H <: AbstractConvexHull = mergehulls!(copy(h), others...)

function merge_hull_lists!(mergedhull::AbstractList, hulltargets::Vector{<:AbstractList}, rev::Bool, orientation::HullOrientation, collinear::Bool, sortedby::Function, targetscollinear::Vector{Bool}, partial::Bool, upper::Bool)  
    # determine starting and stopping points
    f = x -> sortedby(x.data)
    start = rev ? argmax(f, [head(x) for x in hulltargets]) :
                                argmin(f, [head(x) for x in hulltargets])
    stop = start
    if partial
        stop = rev ? argmin(f, [tail(x) for x in hulltargets]) :
                                   argmax(f, [tail(x) for x in hulltargets])
    end

    maxlength = sum(x->length(x), hulltargets)

    # add first point to hull
    empty!(mergedhull)
    pushfirst!(mergedhull, start.data)
    addtarget!(head(mergedhull), start.target)

    maxlength <= 1 && return mergedhull
    
    # prepare orientation test
    betterturn(prevedge,o,a,b) = collinear ? iscloserturn(!orientation,prevedge,o,a,b) : isfurtherturn(!orientation,prevedge,o,a,b)

    # perform jarvis march with search that makes use of the sorted nature of the hulls
    counter = 0
    current = start
    candidates = [head(x) for x in hulltargets]
    prevedge = upper ? UP : DOWN
    while counter == 0 || current !== stop
        if counter > maxlength
            throw(ErrorException("More points were added to the hull than exist in the original hulls to be merged."))
        end
        counter += 1
        for (i, ht) in enumerate(hulltargets)
            candidates[i] = (current.list === ht && (collinear || !targetscollinear[i])) ?      # If the current point belongs to the list being considered, we already know 
                attail(current.next) ? head(ht) : current.next :                                # its candidate point as long as it doesn't contain extraneous collinear points
                jarvissortedsearch(current, prevedge, ht, betterturn)
        end
        next = jarvissearch(current, prevedge, candidates, betterturn)
        if current == next
            throw(ErrorException("Jarvis March failed to progress."))
        end
        # stop adding points when the stopping point has been reached
        if next == stop
            break
        end
        # add the next node to the hull
        push!(mergedhull, next.data)
        addtarget!(tail(mergedhull), next.target)
        prevedge = (next.data[1] - current.data[1], next.data[2] - current.data[2])
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
    hulltargets = [hl.hull for hl in h.subhulls]
    rev = buildinreverse(h.subhulls[1])
    targetscollinear = fill(h.collinear, length(h.subhulls))
    upper = eltype(h.subhulls) <: MutableUpperConvexHull
    partial = upper || eltype(h.subhulls) <: MutableLowerConvexHull
    merge_hull_lists!(mergedhull, hulltargets, rev, h.orientation, h.collinear, h.sortedby, targetscollinear, partial, upper)
    return h
end