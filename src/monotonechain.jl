"""
    lh = lower_monotonechain!(list::PairedLinkedList)

Return the lower convex hull of the points contained in the provided `list` using the [monotone chain algorithm](https://doi.org/10.1016/0020-0190(79)90072-3). 
Each node in the list should contain a two-dimensional point, and the nodes are assumed to be sorted 
(e.g. by lowest "x" value and by lowest "y" in case of ties, though some other sorting methods may produce valid results).
"""
function lower_monotonechain!(pointslist::PairedLinkedList{T}, lower::Union{PairedLinkedList{T},Nothing}=nothing; orientation::HullOrientation = CCW, colinear::Bool = false) where T
    # initialize the convex hull
    if isnothing(lower)
        lower = PairedLinkedList{T}()
        addpartner!(pointslist, lower)
    end
    # exclude or include colinear points on the hull
    wrongturn(args...) = colinear ? !oriented_turn(orientation, args...) : !aligned_turn(orientation, args...)
    # perform monotone chain algorithm
    len = 0
    for node in IteratingListNodes(pointslist; rev=(orientation===CW))
        while len >= 2 && wrongturn(lower.tail.prev.prev.data, lower.tail.prev.data, node.data)
            pop!(lower)
            len -= 1
        end
        push!(lower, node.data)
        if !haspartner(node) # avoid overwriting the partners for points already on the hull
            addpartner!(lower.tail.prev, node)
        end
        len += 1
    end
    return lower
end

"""
    uh = upper_monotonechain!(list::PairedLinkedList)

Return the upper convex hull of the points contained in the provided `list` using the [monotone chain algorithm](https://doi.org/10.1016/0020-0190(79)90072-3). 
Each node in the list should contain a two-dimensional point, and the nodes are assumed to be sorted 
(e.g. by lowest "x" value and by lowest "y" in case of ties, though some other sorting methods may produce valid results).
"""
function upper_monotonechain!(pointslist::PairedLinkedList{T}, upper::Union{PairedLinkedList{T},Nothing}=nothing; orientation::HullOrientation = CCW, colinear::Bool = false) where T
    # initialize the convex hull
    if isnothing(upper)
        upper = PairedLinkedList{T}()
        addpartner!(pointslist, upper)
    end
    # exclude or include colinear points on the hull
    wrongturn(args...) = colinear ? !oriented_turn(orientation, args...) : !aligned_turn(orientation, args...)
    # perform monotone chain algorithm
    len = 0
    for node in IteratingListNodes(pointslist; rev=(orientation===CCW))
        while len >= 2 && wrongturn(upper.tail.prev.prev.data, upper.tail.prev.data, node.data)
            pop!(upper)
            len -= 1
        end
        push!(upper, node.data)
        if !haspartner(node) # avoid overwriting the partners for points already on the hull
            addpartner!(upper.tail.prev, node)
        end
        len += 1
    end
    return upper
end

"""
    h = monotonechain!(list::PairedLinkedList)

Return the convex hull of the points contained in the provided `list` using the [monotone chain algorithm](https://doi.org/10.1016/0020-0190(79)90072-3). 
Each node in the list should contain a two-dimensional point, and the nodes are assumed to be sorted 
(e.g. by lowest "x" value and by lowest "y" in case of ties, though some other sorting methods may produce valid results).
"""
function monotonechain!(pointslist::PairedLinkedList{T}; kwargs...) where T
    # initialize the hull
    hull = PairedLinkedList{T}()
    addpartner!(pointslist, hull)
    # handle the 0- and 1-point cases
    length(pointslist) == 1 && push!(hull, first(pointslist))
    length(pointslist) <= 1 && return hull

    # obtain the lower convex hull
    hull = lower_monotonechain!(pointslist, hull; kwargs...)
    pop!(hull) # remove the last point to avoid duplication (will also be the first point in the upper hull)
    second_from_right = hull.tail.prev

    # obtain the upper convex hull
    hull = upper_monotonechain!(pointslist, hull; kwargs...)
    pop!(hull) # remove the last point to avoid duplication (is also be the first point in the lower hull)

    # re-add partners for the popped nodes
    addpartner!(pointslist.head.next, hull.head.next)
    addpartner!(pointslist.tail.prev, second_from_right.next)
    return hull
end

function lower_monotonechain(points::AbstractVector{T}; orientation::HullOrientation = CCW, colinear::Bool = false, kwargs...) where T
    sortedpoints = sort(points; kwargs...)
    pointslist = PairedLinkedList{T}(sortedpoints...)
    hull = lower_monotonechain!(pointslist; orientation=orientation, colinear=colinear)
    return MutableLowerConvexHull{T}(hull, orientation, colinear)
end

function upper_monotonechain(points::AbstractVector{T}; orientation::HullOrientation = CCW, colinear::Bool = false, kwargs...) where T
    sortedpoints = sort(points; kwargs...)
    pointslist = PairedLinkedList{T}(sortedpoints...)
    hull = upper_monotonechain!(pointslist; orientation=orientation, colinear=colinear)
    return MutableUpperConvexHull{T}(hull, orientation, colinear)
end

function monotonechain(points::AbstractVector{T}; orientation::HullOrientation = CCW, colinear::Bool = false, kwargs...) where T
    sortedpoints = sort(points; kwargs...)
    pointslist = PairedLinkedList{T}(sortedpoints...)
    hull = monotonechain!(pointslist; orientation=orientation, colinear=colinear)
    return MutableConvexHull{T}(hull, orientation, colinear)
end

