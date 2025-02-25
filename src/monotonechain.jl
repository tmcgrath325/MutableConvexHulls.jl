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
    lastdata = start.data
    penultimatedata = lastdata
    for node in ListNodeIterator(start; rev=buildinreverse(h))
        nodedata = node.data
        if !isempty(nodestoadd) && nodestoadd[end] !== node
            while length(nodestoadd) >= 2 # && wrongturn(nodestoadd[end-1].data, nodestoadd[end].data, node.data)
                if coordsareequal(penultimatedata, lastdata) || wrongturn(penultimatedata, lastdata, nodedata)
                    removednode = nodestoadd[end]
                    hastarget(removednode) && deletenode!(removednode.target)
                    pop!(nodestoadd)
                    lastdata = penultimatedata
                    if (length(nodestoadd) > 1)
                        penultimatedata = nodestoadd[end-1].data
                    end
                else
                    break
                end
            end
        end
        push!(nodestoadd, node)
        penultimatedata = lastdata
        lastdata = nodedata
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
    elseif (length(h.points) == 1) || coordsareequal(start.data, stop.data)
        empty!(h.hull)
        insertafter!(newnode(h.hull, head(h.points).data), h.hull.head)
        addtarget!(head(h.hull), head(h.points))
        return h
    end
    # remove extreme points from the hull (allows them to switch ends in the case of certain deletions)
    length(h) > 0 && deletenode!(head(h.hull))    
    length(h) > 0 && deletenode!(tail(h.hull))    
    # exclude or include collinear points on the hull
    wrongturn(o,a,b) = h.collinear ? !isorientedturn(h.orientation,o,a,b) : !isshorterturn(h.orientation,o,a,b)
    # obtain the lower convex hull, and keep track of duplicate coordinates
    lowerhullidxs = Int[]
    lowerhulldups = Bool[]
    hullnode = h.hull.head
    len = 0
    # @show hullnode.data, h.hull
    # println("Forward")
    for (i,node) in enumerate(ListNodeIterator(h.points; rev=buildinreverse(h)))
        o = len == 1 ? hullnode.data : hullnode.prev.data
        a = hullnode.data
        b = node.data
        # println("   consider $(b)")
        if i > 1 && coordsareequal(a,b)
            # println("      ignore equal coords: $(a), $(b)")
            push!(lowerhullidxs, i)
            push!(lowerhulldups, true)
            if hastarget(node)
                deletenode!(node.target)
            end
            continue
        end
        while len >= 2 && wrongturn(o, a, b)
            # println("      remove $(a)")
            while(lowerhulldups[end])
                pop!(lowerhulldups)
                pop!(lowerhullidxs)
            end
            removed = pop!(lowerhullidxs)
            pop!(lowerhulldups)
            # println("         removed idx: $(removed)")
            hullnode = hullnode.prev
            deletenode!(hullnode.next)
            len -= 1
            o = len == 1 ? hullnode.data : hullnode.prev.data
            a = hullnode.data
        end
        if !hastarget(node) # avoid overwriting the targets for points already on the hull
            # println("      add $(b)")
            insertafter!(newnode(h.hull, b), hullnode)
            hullnode = hullnode.next
            addtarget!(hullnode, node)
        else
            if node.target === hullnode.next
                # println("      keep $(hullnode.next.data)")
                hullnode = hullnode.next
            else 
                # println("      move $(node.data) to after $(hullnode.data)")
                movednode = deletenode!(node.target)
                insertafter!(movednode, hullnode)
                addtarget!(movednode, node)
                hullnode = movednode
            end
        end
        push!(lowerhullidxs, i)
        push!(lowerhulldups, false)
        len += 1
        # println("      lower idxs: $(lowerhullidxs)")
    end

    # @show h.hull
    # @show lowerhullidxs
    # @show len
    # println("Backward")

    lidx = length(lowerhullidxs)
    len = 1
    for (i,node) in enumerate(ListNodeIterator(h.points; rev=(h.orientation===CCW)))
        o = len == 1 ? hullnode.data : hullnode.prev.data
        a = hullnode.data
        b = node.data
        # println("   consider $(b)")
        if lidx > 1 && i === h.points.len - lowerhullidxs[lidx] + 1
            lidx -= 1
            if node.target.prev !== hullnode
                # println("      ignore upper hull: $(b)")
                continue
            end
        end 
        if coordsareequal(a,b)
            # println("      ignore equal coords: $(a), $(b)")
            if hastarget(node)
                deletenode!(node.target)
            end
            continue
        end
        while len >= 2 && wrongturn(o, a, b)
            # println("      remove $(a)")
            hullnode = hullnode.prev
            deletenode!(hullnode.next)
            len -= 1
            o = len == 1 ? hullnode.data : hullnode.prev.data
            a = hullnode.data
        end
        if !hastarget(node) # avoid overwriting the targets for points already on the hull
            # println("      add $(b)")
            insertafter!(newnode(h.hull, b), hullnode)
            hullnode = hullnode.next
            addtarget!(hullnode, node)
        else
            # println("      existing hull node $(node.data)")
            if node.target === hullnode.next
                # println("      keep $(hullnode.next.data)")
                hullnode = hullnode.next
            elseif lidx === 1 && i === lowerhullidxs[1]
                # println("      move $(node.data) to after $(hullnode.data)")
                movednode = deletenode!(node.target)
                insertafter!(movednode, hullnode)
                addtarget!(movednode, node)
                hullnode = movednode
            end
        end
        if coordsareequal(hullnode.data, head(h.hull).data)
            deletenode!(hullnode)
            break
        end
        len += 1
    end
    # @show h.hull
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
