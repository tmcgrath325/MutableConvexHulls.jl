firstpoint(h::Union{MutableConvexHull, MutableLowerConvexHull}) = h.orientation === CCW ? h.hull.partner.head.next : h.hull.partner.tail.prev
firstpoint(h::MutableUpperConvexHull) =  h.orientation === CCW ? h.hull.partner.tail.prev : h.hull.partner.head.next

lastpoint(h::Union{MutableConvexHull, MutableLowerConvexHull}) = h.orientation === CCW ? h.hull.partner.tail.prev : h.hull.partner.head.next
lastpoint(h::MutableUpperConvexHull) = h.orientation === CCW ? h.hull.partner.head.next : h.hull.partner.tail.prev

buildinreverse(h::Union{MutableConvexHull, MutableLowerConvexHull}) = h.orientation === CW
buildinreverse(h::MutableUpperConvexHull) = h.orientation === CCW

"""
    lh = lower_monotonechain!(hull [, stop])

Return the lower convex hull of the points contained in the provided `list` using the [monotone chain algorithm](https://doi.org/10.1016/0020-0190(79)90072-3). 
Each node in the list should contain a two-dimensional point, and the nodes are assumed to be sorted 
(e.g. by lowest "x" value and by lowest "y" in case of ties, though some other sorting methods may produce valid results).
"""
function monotonechain!(h::Union{MutableLowerConvexHull, MutableUpperConvexHull},
                        start::PairedListNode{T} = firstpoint(h),
                        stop::PairedListNode{T} = lastpoint(h)) where T
    start.list === stop.list === h.hull.partner || throw(ArgumentError("The start and stop nodes do not belong to the appropriate list."))
    # exclude or include collinear points on the hull
    wrongturn(args...) = h.collinear ? !isorientedturn(h.orientation, args...) : !isshorterturn(h.orientation, args...)
    # perform monotone chain algorithm
    hullnode = haspartner(start) ? start.partner : h.hull.head
    len = length(h.hull.partner) == 0 || start === firstpoint(h) ? 0 : 1
    for node in ListNodeIterator(start; rev=buildinreverse(h))
        if hullnode !== node.partner
            while len >= 2 && wrongturn(hullnode.prev.data, hullnode.data, node.data)
                hullnode = hullnode.prev
                deletenode!(hullnode.next)
                len -= 1
            end
        end
        if !haspartner(node) # avoid overwriting the partners for points already on the hull
            insertnode!(newnode(h.hull, node.data), hullnode)
            hullnode = hullnode.next
            addpartner!(hullnode, node)
        else
            hullnode = node.partner
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
                        start::PairedListNode{T} = firstpoint(h),
                        stop::PairedListNode{T} = lastpoint(h)) where T
    # remove extreme points from the hull (allows them to switch ends in the case of certain deletions)
    length(h) > 0 && deletenode!(head(h.hull))    
    length(h) > 0 && deletenode!(tail(h.hull))    
    # exclude or include collinear points on the hull
    wrongturn(args...) = h.collinear ? !isorientedturn(h.orientation, args...) : !isshorterturn(h.orientation, args...)
    # obtain the lower convex hull
    onlowerhull = fill(false, h.hull.partner.len)
    hullnode = h.hull.head
    len = 0
    for (i,node) in enumerate(ListNodeIterator(h.hull.partner; rev=buildinreverse(h)))
        while len >= 2 && wrongturn(hullnode.prev.data, hullnode.data, node.data)
            removedindex = findlast(onlowerhull)
            onlowerhull[removedindex] = false
            hullnode = hullnode.prev
            deletenode!(hullnode.next)
            len -= 1
        end
        if !haspartner(node) # avoid overwriting the partners for points already on the hull
            insertnode!(newnode(h.hull, node.data), hullnode)
            hullnode = hullnode.next
            addpartner!(hullnode, node)
        else
            if node.partner != hullnode.next
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
    for (node, islower) in zip(ListNodeIterator(h.hull.partner; rev=(h.orientation===CCW)), onlowerhull)
        islower && continue
        while len >= 2 && wrongturn(hullnode.prev.data, hullnode.data, node.data)
            hullnode = hullnode.prev
            deletenode!(hullnode.next)
            len -= 1
        end
        if !haspartner(node) # avoid overwriting the partners for points already on the hull
            addnode = newnode(h.hull, node.data)
            insertnode!(addnode, hullnode)
            hullnode = hullnode.next
            addpartner!(hullnode, node)
        else
            hullnode = hullnode.next
        end
        len += 1
    end

    return h
end

function lower_monotonechain(points::AbstractVector{T}; orientation::HullOrientation = CCW, collinear::Bool = false, sortedby::Function = identity, presorted::Bool = false) where T
    sortedpoints = presorted ? points : sort(points; by = sortedby)
    pointslist = PairedLinkedList{T}(sortedpoints...)
    hull = PairedLinkedList{T}()
    addpartner!(hull, pointslist)
    h = MutableLowerConvexHull{T, typeof(sortedby)}(hull, orientation, collinear, sortedby, true)
    monotonechain!(h)
    return h
end
lower_monotonechain(points::NTuple{N,T}; kwargs...) where {N,T} = lower_monotonechain([points...]; kwargs...)

function upper_monotonechain(points::AbstractVector{T}; orientation::HullOrientation = CCW, collinear::Bool = false, sortedby::Function = identity, presorted::Bool = false) where T
    sortedpoints = presorted ? points : sort(points; by = sortedby)
    pointslist = PairedLinkedList{T}(sortedpoints...)
    hull = PairedLinkedList{T}()
    addpartner!(hull, pointslist)
    h = MutableUpperConvexHull{T, typeof(sortedby)}(hull, orientation, collinear, sortedby, true)
    monotonechain!(h)
    return h
end
upper_monotonechain(points::NTuple{N,T}; kwargs...) where {N,T} = upper_monotonechain([points...]; kwargs...)

function monotonechain(points::AbstractVector{T}; orientation::HullOrientation = CCW, collinear::Bool = false, sortedby::Function = identity, presorted::Bool = false) where T
    sortedpoints = presorted ? points : sort(points; by = sortedby)
    pointslist = PairedLinkedList{T}(sortedpoints...)
    hull = PairedLinkedList{T}()
    addpartner!(hull, pointslist)
    h = MutableConvexHull{T, typeof(sortedby)}(hull, orientation, collinear, sortedby, true)
    monotonechain!(h)
    return h
end
monotonechain(points::NTuple{N,T}; kwargs...) where {N,T} = monotonechain([points...]; kwargs...)

