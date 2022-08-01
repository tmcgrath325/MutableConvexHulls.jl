# cross product of 2D vectors OA and OB
cross(o, a, b) = (a[1] - o[1]) * (b[2] - o[2]) - (a[2] - o[2]) * (b[1] - o[1])

"""
lh = lower_grahamscan(list::PairedLinkedList)

Return the lower convex hull of the points contained in the provided `list`. Each node in the list should contain
a two-dimensional point, and the nodes are assumed to be sorted (e.g. by lowest "x" value and by lowest "y" in case of ties).
"""
function lower_grahamscan(pointsList::PairedLinkedList{T}) where T
    # initialize the lower convex hull
    lower = PairedLinkedList{T}(pointsList)
    pointsList.partner = lower
    lower.head.partner = pointslist.head
    pointslist.head.partner = lower.head
    lower.tail.partner = pointslist.tail
    pointslist.tail.partner = lower.tail
    # perform monotone chain algorithm
    pointsNode = pointsList.head
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

"""
lh = lower_grahamscan(list::PairedLinkedList)

Return the lower convex hull of the points contained in the provided `list`. Each node in the list should contain
a two-dimensional point, and the nodes are assumed to be sorted (e.g. by lowest "x" value and by lowest "y" in case of ties).
"""
function upper_grahamscan(pointsList::PairedLinkedList{T}) where T
    # initialize the lower convex hull
    upper = PairedLinkedList{T}(pointsList)
    pointsList.partner = upper
    upper.head.partner = pointslist.head
    pointslist.head.partner = upper.head
    upper.tail.partner = pointslist.tail
    pointslist.tail.partner = upper.tail
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

"""
h = grahamscan(list::PairedLinkedList)

Return the convex hull of the points contained in the provided `list`. Each node in the list should contain
a two-dimensional point, and the nodes are assumed to be sorted (e.g. by lowest "x" value and by lowest "y" in case of ties).
"""
function grahamscan(pointsList::PairedLinkedList{T}) where T
    # obtain upper and lower hulls
    upper = upper_hull(pointsList)
    lower = lower_hull(pointsList)
    # stitch the two together, without repeating the far left and far right points
    pop!(upper)
    popfirst!(upper)
    append!(lower, upper)
    return hull
end

