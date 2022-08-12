function jarvis_binary_search(query::AbstractListNode, hull::AbstractConvexHull)
    # if the query explicitly belongs to the hull, we already know the appropriate next point
    query.list === hullpointslist && return (PairedLinkedLists.at_tail(query.next) ? pointslist.head.next : query.next)
    return jarvis_binary_search(query, hull.hull, hull.collinear, hull.orientation)
end

function jarvis_binary_search(query::AbstractListNode, hullpointslist::AbstractLinkedList, collinear::Bool = false, orientation::HullOrientation = CCW)  
    betterturn(args...) = collinear ? closer_turn(!orientation, args...) : further_turn(!orientation, args...)
    return jarvis_binary_search(query, hullpointslist, betterturn)
end

function jarvis_binary_search(query::AbstractListNode, pointslist::AbstractLinkedList, betterturn::Function)
    pointslist.len == 0 && throw(ArgumentError("The list of points must not be empty."))
    pointslist.len == 1 && return pointslist.head.next
    pointslist.head.next == pointslist.tail.prev && throw(ArgumentError("All points in the list are duplicates."))

    # initialize interval bounds
    left = 1
    right = pointslist.len
    middle = Int(floor(pointslist.len/2))
    target = getnode(pointslist, middle)

    # recursively shrink interval until a tangent is found
    return jarvis_binary_search_interval(left, right, middle, target, query, pointslist, betterturn)
end

# Perform binary search in the interval specified by left and right. The middle index and target node are already determined for the interval
function jarvis_binary_search_interval(left::Int, right::Int, middle::Int, target::AbstractListNode, query::AbstractListNode, pointslist::AbstractLinkedList, betterturn::Function)
    # wrap around the list to get prev/next values for the ends
    prev = middle == 1 ? pointslist.tail.prev : target.prev
    next = middle == pointslist.len ? pointlist.head.next : target.next

    # if the query is a duplicate of a hull point, 
    (query.data == target.data) && return next
    (query.data == prev.data) && return target
    (query.data == next.data) && return next.next

    # check if the adjacent nodes represent "better" options than the target one
    prev_better = betterturn(query.data, target.data, next.data)
    next_better = betterturn(query.data, target.data, prev.data)
    if prev_better && next_better # the target is the best choice
        return target
    elseif prev_better # explore the left interval
        (right - left) == 0 && return prev
        right = middle
        middle = Int(floor((right - left)/2) + left)
        idx_change = right - middle
        for i=1:idx_change
            target = target.prev
        end
    else # explore the right interval
        (right - left) == 0 && return next
        left = middle
        middle = Int(floor((right - left)/2) + left)
        idx_change = middle - left
        for i=1:idx_change
            target = target.next
        end
    end

    return jarvis_binary_search_interval(left, right, middle, target, query, pointslist, betterturn)
end

function tangent(query::AbstractListNode, hullpointslist::AbstractLinkedList, betterturn::Function)

end