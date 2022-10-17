""" 
    node = firstpoint(h)

Get the first point on `h.points` that should be considered for the monotone chain algorithm.
"""
firstpoint(h::Union{MutableConvexHull, MutableLowerConvexHull}) = h.orientation === CCW ? h.points.head.next : h.points.tail.prev
firstpoint(h::MutableUpperConvexHull) =  h.orientation === CCW ? h.points.tail.prev : h.points.head.next

""" 
    node = lastpoint(h)

Get the last point on `h.points` that should be considered for the monotone chain algorithm.
"""
lastpoint(h::Union{MutableConvexHull, MutableLowerConvexHull}) = h.orientation === CCW ? h.points.tail.prev : h.points.head.next
lastpoint(h::MutableUpperConvexHull) = h.orientation === CCW ? h.points.head.next : h.points.tail.prev

""" 
    buildinreverse(h) -> Bool

Return `true` if the monotonechain algorithm should start at the last point in `h.points`, and `false` otherwise.
"""
buildinreverse(h::Union{MutableConvexHull, MutableLowerConvexHull}) = h.orientation === CW
buildinreverse(h::MutableUpperConvexHull) = h.orientation === CCW

"""
    monotonechain!(hull [, start, stop])

Determine the convex hull of the points contained in the provided `hull.points` using the [monotone chain algorithm](https://doi.org/10.1016/0020-0190(79)90072-3). 
Each node in the list should contain a two-dimensional point, and the nodes are assumed to be sorted 
(e.g. by lowest "x" value and by lowest "y" in case of ties, though some other sorting methods may produce valid results).

`start` and `stop` should be nodes contained in `hull.points`. If `start` and/or `stop` are provided, the convex hull will only be updated on the inclusive interval
between `start` and `stop`. This allows efficient removal of points from a convex hull.
"""
function monotonechain!(h::Union{MutableLowerConvexHull, MutableUpperConvexHull},
                        start::PointNode{T} = firstpoint(h),
                        stop::PointNode{T} = lastpoint(h)) where T
    if isempty(h.points)
        empty!(h.hull)
        return h
    end
    start.list === stop.list === h.points || throw(ArgumentError("The start and stop nodes do not belong to the appropriate list."))
    # exclude or include collinear points on the hull
    wrongturn(o,a,b) = h.collinear ? !isorientedturn(h.orientation,o,a,b) : !isshorterturn(h.orientation,o,a,b)
    # get a list of point nodes to be added to the hull
    nodestoadd = PointNode{T}[]
    start === firstpoint(h) && push!(nodestoadd, start)
    for node in ListNodeIterator(start; rev=buildinreverse(h))
        if !isempty(nodestoadd) && nodestoadd[end] !== node
            while length(nodestoadd) >= 2 # && wrongturn(nodestoadd[end-1].data, nodestoadd[end].data, node.data)
                if coordsareequal(nodestoadd[end-1].data, nodestoadd[end].data) || wrongturn(nodestoadd[end-1].data, nodestoadd[end].data, node.data)
                    removednode = nodestoadd[end]
                    hastarget(removednode) && deletenode!(removednode.target)
                    pop!(nodestoadd)
                else
                    break
                end
            end
        end
        push!(nodestoadd, node)
        if node === stop
            if length(nodestoadd) >= 2 && coordsareequal(nodestoadd[end-1].data, stop.data)
                removednode = nodestoadd[end]
                hastarget(removednode) && deletenode!(removednode.target)
                pop!(nodestoadd)
            end
            break
        end
    end
    # add the appropriate point nodes to the hull
    hullnode = hastarget(start) ? start.target : h.hull.head
    for node in nodestoadd
        if !hastarget(node) # avoid overwriting the targets for points already on the hull
            insertafter!(newnode(h.hull, node.data), hullnode)
            hullnode = hullnode.next
            addtarget!(hullnode, node)
        else
            hullnode = node.target
        end
    end
    return h
end

# TO DO: implement use of start and stop arguments, and avoid unecessary node constructor calls
function monotonechain!(h::MutableConvexHull{T},
                        start::PointNode{T} = firstpoint(h),
                        stop::PointNode{T} = lastpoint(h)) where T
    if isempty(h.points)
        empty!(h.hull)
        return h
    end
    # remove extreme points from the hull (allows them to switch ends in the case of certain deletions)
    length(h) > 0 && deletenode!(head(h.hull))    
    length(h) > 0 && deletenode!(tail(h.hull))    
    # exclude or include collinear points on the hull
    wrongturn(o,a,b) = h.collinear ? !isorientedturn(h.orientation,o,a,b) : !isshorterturn(h.orientation,o,a,b)
    # obtain the lower convex hull
    lowerhullidxs = Int[]
    hullnode = h.hull.head
    len = 0
    for (i,node) in enumerate(ListNodeIterator(h.points; rev=buildinreverse(h)))
        if len == 0 || !coordsareequal(hullnode.data, node.data)
            while len >= 2 && wrongturn(hullnode.prev.data, hullnode.data, node.data)
                pop!(lowerhullidxs)
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
            i !== 1 && push!(lowerhullidxs, i)
            len += 1
        end
    end

    lidx = length(lowerhullidxs)
    len = 1
    for (i,node) in enumerate(ListNodeIterator(h.points; rev=(h.orientation===CCW)))
        if lidx > 0 && i === h.points.len - lowerhullidxs[lidx] + 1
            lidx -= 1
            continue
        end 
        if !coordsareequal(hullnode.data, node.data)
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
    end

    if length(h) > 1 && coordsareequal(head(h.hull).data, tail(h.hull).data)
        deletenode!(tail(h.hull))
    end

    return h
end

"""
    lh = lower_monotonechain(points [; orientation, collinear, sortedby])

Return the lower convex hull generated from the provided `points`.

`orientation` specifies whether the points along the convex hull are ordered clockwise `CW`, or counterclockwise `CCW`, and defaults to `CCW`.

`collinear` specifies whether collinear points are allowed along the surface of the convex hull, and defaults to `false`.

`sortedby` specifies a function to apply to points prior to sorting, and defaults to `identity` (resulting in default sorting behavior).
"""
function lower_monotonechain(points::AbstractVector{T}; orientation::HullOrientation = CCW, collinear::Bool = false, sortedby::Function = identity) where T
    pointslist = PointList{T}(;sortedby=sortedby)
    hull = HullList{T,typeof(sortedby)}()
    addtarget!(hull, pointslist)
    for p in points
        push!(pointslist,p)
    end
    h = MutableLowerConvexHull{T, typeof(sortedby)}(hull, pointslist, orientation, collinear, sortedby)
    monotonechain!(h)
    return h
end
lower_monotonechain(points::Matrix; kwargs...) = lower_monotonechain([(points[i,:]...,) for i=1:size(points,1)])

"""
    uh = upper_monotonechain(points [; orientation, collinear, sortedby])

Return the upper convex hull generated from the provided `points`.

`orientation` specifies whether the points along the convex hull are ordered clockwise `CW`, or counterclockwise `CCW`, and defaults to `CCW`.

`collinear` specifies whether collinear points are allowed along the surface of the convex hull, and defaults to `false`.

`sortedby` specifies a function to apply to points prior to sorting, and defaults to `identity` (resulting in default sorting behavior).
"""
function upper_monotonechain(points::AbstractVector{T}; orientation::HullOrientation = CCW, collinear::Bool = false, sortedby::Function = identity) where T
    pointslist = PointList{T}(;sortedby=sortedby)
    hull = HullList{T,typeof(sortedby)}()
    addtarget!(hull, pointslist)
    for p in points
        push!(pointslist,p)
    end
    h = MutableUpperConvexHull{T, typeof(sortedby)}(hull, pointslist, orientation, collinear, sortedby)
    monotonechain!(h)
    return h
end
upper_monotonechain(points::Matrix; kwargs...) = upper_monotonechain([(points[i,:]...,) for i=1:size(points,1)])

"""
    h = monotonechain(points [; orientation, collinear, sortedby])

Return the convex hull generated from the provided `points`.

`orientation` specifies whether the points along the convex hull are ordered clockwise `CW`, or counterclockwise `CCW`, and defaults to `CCW`.

`collinear` specifies whether collinear points are allowed along the surface of the convex hull, and defaults to `false`.

`sortedby` specifies a function to apply to points prior to sorting, and defaults to `identity` (resulting in default sorting behavior).
"""
function monotonechain(points::AbstractVector{T}; orientation::HullOrientation = CCW, collinear::Bool = false, sortedby::Function = identity) where T
    pointslist = PointList{T}(;sortedby=sortedby)
    hull = HullList{T,typeof(sortedby)}()
    addtarget!(hull, pointslist)
    for p in points
        push!(pointslist,p)
    end
    h = MutableConvexHull{T, typeof(sortedby)}(hull, pointslist, orientation, collinear, sortedby)
    monotonechain!(h)
    return h
end
monotonechain(points::Matrix; kwargs...) = monotonechain([(points[i,:]...,) for i=1:size(points,1)])
