function jarvissortedsearch(query::AbstractListNode, prevedge, pointslist::AbstractLinkedList, betterturn::Function)
    pointslist.len == 0 && throw(ArgumentError("The list of points must not be empty."))
    pointslist.len == 1 && return pointslist.head.next
    pointslist.head.next == pointslist.tail.prev && throw(ArgumentError("All points in the list are duplicates."))

    prev_worse = betterturn(prevedge, query.data, pointslist.tail.prev.data, pointslist.head.next.data)
    for target in IteratingListNodes(pointslist)
        if query.data == target.next.data
            continue
        end
        next_worse = !betterturn(prevedge, query.data, target.data, target.next.data)
        if prev_worse && next_worse
            return target
        end
        prev_worse = !next_worse
    end
    return pointslist.tail.prev
end

# function jarvisbinarysearch(query::AbstractListNode, prevedge, pointslist::AbstractLinkedList, betterturn::Function)
#     pointslist.len == 0 && throw(ArgumentError("The list of points must not be empty."))
#     pointslist.len == 1 && return pointslist.head.next
#     pointslist.head.next == pointslist.tail.prev && throw(ArgumentError("All points in the list are duplicates."))

#     # initialize interval bounds
#     left = 1
#     right = pointslist.len
#     middle = Int(floor(pointslist.len/2))
#     target = getnode(pointslist, middle)

#     # recursively shrink interval until a tangent is found
#     return jarvisbinarysearchinterval(left, right, middle, target, query, prevedge, pointslist, betterturn)
# end

# # Perform binary search in the interval specified by left and right. The middle index and target node are already determined for the interval
# function jarvisbinarysearchinterval(left::Int, right::Int, middle::Int, target::AbstractListNode, query::AbstractListNode, prevedge, pointslist::AbstractLinkedList, betterturn::Function)
#     left == middle == right && return target

#     # wrap around the list to get prev/next values for the ends
#     prev = middle == 1 ? pointslist.tail.prev : target.prev
#     next = middle == pointslist.len ? pointslist.head.next : target.next

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

function mergehulls(h::H, others::H...; orientation::HullOrientation=h.orientation, collinear::Bool=h.collinear, by::Function=identity) where {T,H<:AbstractConvexHull{T}}
    # filter out empty hulls
    hulls = filter(x->length(x.hull)>0,[h, others...])
    length(hulls) == 0 && return copy(h)
    length(hulls) == 1 && return copy(only(hulls))

    # prepare sorting function and orientation test
    f = x -> by(x.data)
    betterturn(args...) = collinear ? iscloserturn(!orientation, args...) : isfurtherturn(!orientation, args...)

    # initialize new list of points and new hull
    newpoints = PairedLinkedList{T}()
    newhull = PairedLinkedList{T}()
    addpartner!(newpoints, newhull)

    # Set up copies of the hulls that point to the new points list
    hulltargets = [TargetedLinkedList(newpoints) for i=1:length(hulls)]
    for (i, originalhull) in enumerate(hulls)
        for data in originalhull.hull
            push!(hulltargets[i], data)
        end
    end

    # add points from all hulls into the new points list. 
    # Assuming the list of points for each hull are sorted the same way, this will preserve the sorting in the new hull.
    # Otherwise, there is no guarantee of sorting for the list of points for the new hull. 
    # The hull itself will still be sorted as specified by `orientation` and `collinear`
    remaininghulls = collect(1:length(hulls))
    reversehulls = [x.orientation != orientation for x in hulls]
    pointnodes = [reversehulls[i] ? x.hull.partner.tail.prev : x.hull.partner.head.next  for (i,x) in enumerate(hulls)]
    while !isempty(remaininghulls)
        hullidx = argmin(x->f(pointnodes[x]), remaininghulls)
        nextnode = pointnodes[hullidx]
        push!(newpoints, nextnode.data)
        if haspartner(nextnode)
            targetingnode = getfirst(x->x.data == nextnode.data, IteratingListNodes(hulltargets[hullidx]))
            addpartner!(targetingnode, newpoints.tail.prev)
        end
        pointnodes[hullidx] = reversehulls[hullidx] ? pointnodes[hullidx].prev : pointnodes[hullidx].next
        if reversehulls[hullidx] ? PairedLinkedLists.at_head(pointnodes[hullidx]) : PairedLinkedLists.at_tail(pointnodes[hullidx])
            deleteidx = findfirst(x->x==hullidx, remaininghulls)
            deleteat!(remaininghulls, deleteidx)
        end
    end
    
    # determine starting and stopping points
    start = argmin(f, map(x->argmin(f, IteratingListNodes(x)), hulltargets))
    start = orientation === CCW ? argmin(f, map(x->argmin(f, IteratingListNodes(x)), hulltargets)) : argmax(f, map(x->argmax(f, IteratingListNodes(x)), hulltargets))
    stop = start
    if H <: Union{MutableUpperConvexHull, MutableLowerConvexHull}
        stop = orientation === CCW ? argmax(f, map(x->argmax(f, IteratingListNodes(x)), hulltargets)) : argmin(f, map(x->argmin(f, IteratingListNodes(x)), hulltargets))
    end

    # add first point to hull
    push!(newhull, start.data)
    addpartner!(newhull.tail.prev, start.partner)
    
    # perform jarvis march with search that makes use of the sorted nature of the hulls
    counter = 0
    maxlength = sum(x->length(x.hull), hulls)
    current = newhull.head.next
    candidates = [x.head.next for x in hulltargets]
    prevedge = H <: MutableUpperConvexHull ? UP : DOWN
    while counter == 0 || current !== stop.partner
        if counter > maxlength
            throw(ErrorException("More points were added to the hull than exist in the original hulls to be merged."))
        end
        counter += 1
        for (i, ht) in enumerate(hulltargets)
            candidates[i] = (current.list === ht && (collinear || !hulls[i].collinear)) ?       # if the current point belongs to the list being considered, we already know 
                (reversehulls[i] ?                                                              # its candidate point as long as it doesn't contain extraneous collinear points
                    (PairedLinkedLists.at_head(current.prev) ? ht.tail.prev : current.prev) : 
                    (PairedLinkedLists.at_tail(current.next) ? ht.head.next : current.next)) : 
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
        push!(newhull, next.data)
        addpartner!(newhull.tail.prev, next.partner)
        prevedge = next.data .- current.data
        current = next
    end

    # if these are only partial convex hulls (i.e. upper or lower), add the stopping point at the end of the hull
    if H <: Union{MutableUpperConvexHull, MutableLowerConvexHull}
        push!(newhull, stop.data)
        addpartner!(newhull.tail.prev, stop.partner)
    end

    return H(newhull, orientation, collinear)
end