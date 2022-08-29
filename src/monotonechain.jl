firstpoint(h::Union{MutableConvexHull, MutableLowerConvexHull}) = h.orientation === CCW ? h.points.head.next : h.points.tail.prev
firstpoint(h::MutableUpperConvexHull) =  h.orientation === CCW ? h.points.tail.prev : h.points.head.next

lastpoint(h::Union{MutableConvexHull, MutableLowerConvexHull}) = h.orientation === CCW ? h.points.tail.prev : h.points.head.next
lastpoint(h::MutableUpperConvexHull) = h.orientation === CCW ? h.points.head.next : h.points.tail.prev

buildinreverse(h::Union{MutableConvexHull, MutableLowerConvexHull}) = h.orientation === CW
buildinreverse(h::MutableUpperConvexHull) = h.orientation === CCW

"""
    lh = lower_monotonechain!(hull [, stop])

Return the lower convex hull of the points contained in the provided `list` using the [monotone chain algorithm](https://doi.org/10.1016/0020-0190(79)90072-3). 
Each node in the list should contain a two-dimensional point, and the nodes are assumed to be sorted 
(e.g. by lowest "x" value and by lowest "y" in case of ties, though some other sorting methods may produce valid results).
"""
function monotonechain!(h::Union{MutableLowerConvexHull, MutableUpperConvexHull},
                        start::PointNode{T} = firstpoint(h),
                        stop::PointNode{T} = lastpoint(h)) where T
    start.list === stop.list === h.points || throw(ArgumentError("The start and stop nodes do not belong to the appropriate list."))
    # exclude or include collinear points on the hull
    wrongturn(args...) = h.collinear ? !isorientedturn(h.orientation, args...) : !isshorterturn(h.orientation, args...)
    # perform monotone chain algorithm
    hullnode = hastarget(start) ? start.target : h.hull.head
    len = length(h.points) == 0 || start === firstpoint(h) ? 0 : 1
    for node in ListNodeIterator(start; rev=buildinreverse(h))
        if hullnode !== node.target
            while len >= 2 && wrongturn(hullnode.prev.data, hullnode.data, node.data)
                hullnode = hullnode.prev
                deletenode!(hullnode.next)
                len -= 1
            end
        end
        if !hastarget(node) # avoid overwriting the targets for points already on the hull
            insertafter!(newnode(h.hull, node.data), hullnode)
            hullnode = hullnode.next
            addtarget!(hullnode, node)
        else
            hullnode = node.target
        end
        node === stop && break
        len += 1
    end
    return h
end

"""
    h = monotonechain!(list::PairedLinkedList)

Return the convex hull of the points contained in the provided `list` using the [monotone chain algorithm](https://doi.org/10.1016/0020-0190(79)90072-3). 
Each node in the list should contain a two-dimensional point, and the nodes are assumed to be sorted 
(e.g. by lowest "x" value and by lowest "y" in case of ties, though some other sorting methods may produce valid results).
"""
function monotonechain!(h::MutableConvexHull{T},
                        start::PointNode{T} = firstpoint(h),
                        stop::PointNode{T} = lastpoint(h)) where T
    # remove extreme points from the hull (allows them to switch ends in the case of certain deletions)
    length(h) > 0 && deletenode!(head(h.hull))    
    length(h) > 0 && deletenode!(tail(h.hull))    
    # exclude or include collinear points on the hull
    wrongturn(args...) = h.collinear ? !isorientedturn(h.orientation, args...) : !isshorterturn(h.orientation, args...)
    # obtain the lower convex hull
    onlowerhull = fill(false, h.points.len)
    hullnode = h.hull.head
    len = 0
    for (i,node) in enumerate(ListNodeIterator(h.points; rev=buildinreverse(h)))
        while len >= 2 && wrongturn(hullnode.prev.data, hullnode.data, node.data)
            removedindex = findlast(onlowerhull)
            onlowerhull[removedindex] = false
            hullnode = hullnode.prev
            deletenode!(hullnode.next)
            len -= 1
        end
        if !hastarget(node) # avoid overwriting the targets for points already on the hull
            insertafter!(newnode(h.hull, node.data), hullnode)
            hullnode = hullnode.next
            addtarget!(hullnode, node)
        else
            if node.target != hullnode.next
                continue
            else
                hullnode = hullnode.next
            end
        end
        onlowerhull[i] = i !== 1
        len += 1
    end

    reverse!(onlowerhull)
    # obtain the upper convex hull
    len = 1
    for (node, islower) in zip(ListNodeIterator(h.points; rev=(h.orientation===CCW)), onlowerhull)
        islower && continue
        while len >= 2 && wrongturn(hullnode.prev.data, hullnode.data, node.data)
            hullnode = hullnode.prev
            deletenode!(hullnode.next)
            len -= 1
        end
        if !hastarget(node) # avoid overwriting the targets for points already on the hull
            addnode = newnode(h.hull, node.data)
            insertafter!(addnode, hullnode)
            hullnode = hullnode.next
            addtarget!(hullnode, node)
        else
            hullnode = hullnode.next
        end
        len += 1
    end

    return h
end

function lower_monotonechain(points::AbstractVector{T}; orientation::HullOrientation = CCW, collinear::Bool = false, sortedby::Function = identity) where T
    pointslist = PointList{T}(;sortedby=sortedby)
    hull = HullList{T,typeof(sortedby)}()
    addtarget!(hull, pointslist)
    # push!(pointslist, points...)
    for p in points
        push!(pointslist,p)
    end
    h = MutableLowerConvexHull{T, typeof(sortedby)}(hull, pointslist, orientation, collinear, sortedby)
    monotonechain!(h)
    return h
end
lower_monotonechain(points::NTuple{N,T}; kwargs...) where {N,T} = lower_monotonechain([points...]; kwargs...)

function upper_monotonechain(points::AbstractVector{T}; orientation::HullOrientation = CCW, collinear::Bool = false, sortedby::Function = identity) where T
    pointslist = PointList{T}(;sortedby=sortedby)
    hull = HullList{T,typeof(sortedby)}()
    addtarget!(hull, pointslist)
    # push!(pointslist, points...)
    for p in points
        push!(pointslist,p)
    end
    h = MutableUpperConvexHull{T, typeof(sortedby)}(hull, pointslist, orientation, collinear, sortedby)
    monotonechain!(h)
    return h
end
upper_monotonechain(points::NTuple{N,T}; kwargs...) where {N,T} = upper_monotonechain([points...]; kwargs...)

function monotonechain(points::AbstractVector{T}; orientation::HullOrientation = CCW, collinear::Bool = false, sortedby::Function = identity) where T
    pointslist = PointList{T}(;sortedby=sortedby)
    hull = HullList{T,typeof(sortedby)}()
    addtarget!(hull, pointslist)
    # push!(pointslist, points...)
    h = MutableConvexHull{T, typeof(sortedby)}(hull, pointslist, orientation, collinear, sortedby)
    monotonechain!(h)
    return h
end
monotonechain(points::NTuple{N,T}; kwargs...) where {N,T} = monotonechain([points...]; kwargs...)

