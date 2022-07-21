# cross product of 2D vectors OA and OB
cross(o, a, b) = (a[1] - o[1]) * (b[2] - o[2]) - (a[2] - o[2]) * (b[1] - o[1])

function lower_hull(pointsList::DoublyLinkedList{T}) where T
    # initialize the lower convext hull
    lower = DoublyLinkedList{T}()
    lower.node.pair = pointsList.node.pair
    pointsNode = pointsList.node
    # perform monotone chain algorithm
    for i=1:pointsList.len
        pointsNode = pointsNode.next
        while length(lower) >= 2 && cross(lower.node.prev.prev.data, lower.node.prev.data, pointsNode.data) <= 0
            pop!(lower)
        end
        push!(lower, pointsNode.data)
        lower.node.prev.pair = pointsNode
        pointsNode.pair = lower.node.prev
    end
    return lower
end

function upper_hull(pointsList::DoublyLinkedList{T}) where T
    # initialize the lower convext hull
    upper = DoublyLinkedList{T}()
    upper.node.pair = pointsList.node.pair
    pointsNode = pointsList.node
    # perform monotone chain algorithm
    for i=1:pointsList.len
        pointsNode = pointsNode.prev
        while length(upper) >= 2 && cross(upper.node.prev.prev.data, upper.node.prev.data, pointsNode.data) <= 0
            pop!(upper)
        end
        push!(upper, pointsNode.data)
        upper.node.prev.pair = pointsNode
        pointsNode.pair = upper.node.prev
    end
    return upper
end

function convex_hull(pointsList::DoublyLinkedList{T}) where T
    # obtain upper and lower hulls
    hull = lower_hull(pointsList)
    upper = upper_hull(pointsList)
    # stitch the two together, without repeating the far left and far right points
    hull.node.prev.next = upper.node.next.next
    hull.node.prev = upper.node.prev.prev
    hull.node.prev.next = hull.node  
    hull.len = hull.len + upper.len - 2
    return hull
end

