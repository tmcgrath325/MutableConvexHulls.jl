"""
    lh = lower_monotonechain!(hull [, stop])

Return the lower convex hull of the points contained in the provided `list` using the [monotone chain algorithm](https://doi.org/10.1016/0020-0190(79)90072-3). 
Each node in the list should contain a two-dimensional point, and the nodes are assumed to be sorted 
(e.g. by lowest "x" value and by lowest "y" in case of ties, though some other sorting methods may produce valid results).
"""
function monotonechain!(h::MutableLowerConvexHull, stop::PairedListNode{T} = h.orientation === CCW ? h.hull.partner.tail.prev : h.hull.partner.head.next) where T
    # exclude or include collinear points on the hull
    wrongturn(args...) = h.collinear ? !isorientedturn(h.orientation, args...) : !isshorterturn(h.orientation, args...)
    # perform monotone chain algorithm
    len = 0
    lowernode = h.hull.head
    for node in ListNodeIterator(h.hull.partner; rev=(h.orientation===CW))
        while len >= 2 && wrongturn(lowernode.prev.data, lowernode.data, node.data)
            lowernode = lowernode.prev
            deletenode!(lowernode.next)
            len -= 1
        end
        if !haspartner(node) # avoid overwriting the partners for points already on the hull
            insertnode!(newnode(h.hull, node.data), lowernode)
            lowernode = lowernode.next
            addpartner!(lowernode, node)
        else
            lowernode = node.partner
        end
        node === stop && break
        len += 1
    end
    return h
end

"""
    uh = upper_monotonechain!(hull [, stop])

Return the upper convex hull of the points contained in the provided `list` using the [monotone chain algorithm](https://doi.org/10.1016/0020-0190(79)90072-3). 
Each node in the list should contain a two-dimensional point, and the nodes are assumed to be sorted 
(e.g. by lowest "x" value and by lowest "y" in case of ties, though some other sorting methods may produce valid results).
"""
function monotonechain!(h::MutableUpperConvexHull, stop::PairedListNode{T} = h.orientation === CCW ? h.hull.partner.head.next : h.hull.partner.tail.prev) where T
    # exclude or include collinear points on the hull
    wrongturn(args...) = h.collinear ? !isorientedturn(h.orientation, args...) : !isshorterturn(h.orientation, args...)
    # perform monotone chain algorithm
    len = 0
    uppernode = haspartner(h.hull.partner.tail.prev) ? h.hull.partner.tail.prev.partner : h.hull.head
    for node in ListNodeIterator(h.hull.partner; rev=(h.orientation===CCW))
        if uppernode !== node.partner
            while len >= 2 && wrongturn(uppernode.prev.data, uppernode.data, node.data)
                uppernode = uppernode.prev
                deletenode!(uppernode.next)
                len -= 1
            end
        end
        if !haspartner(node) # avoid overwriting the partners for points already on the hull
            addnode = newnode(h.hull, node.data)
            insertnode!(addnode, uppernode)
            uppernode = uppernode.next
            addpartner!(uppernode, node)
        else
            uppernode = node.partner
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
function monotonechain!(h::MutableConvexHull{T}) where T
    # exclude or include collinear points on the hull
    wrongturn(args...) = h.collinear ? !isorientedturn(h.orientation, args...) : !isshorterturn(h.orientation, args...)

    # obtain the lower convex hull
    onlowerhull = fill(false, h.hull.partner.len)
    hullnode = h.hull.head
    len = 0
    for (i,node) in enumerate(ListNodeIterator(h.hull.partner; rev=(h.orientation===CW)))
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
    h = MutableLowerConvexHull{T}(hull, orientation, collinear, sortedby)
    monotonechain!(h)
    return h
end

function upper_monotonechain(points::AbstractVector{T}; orientation::HullOrientation = CCW, collinear::Bool = false, sortedby::Function = identity, presorted::Bool = false) where T
    sortedpoints = presorted ? points : sort(points; by = sortedby)
    pointslist = PairedLinkedList{T}(sortedpoints...)
    hull = PairedLinkedList{T}()
    addpartner!(hull, pointslist)
    h = MutableUpperConvexHull{T}(hull, orientation, collinear, sortedby)
    monotonechain!(h)
    return h
end

function monotonechain(points::AbstractVector{T}; orientation::HullOrientation = CCW, collinear::Bool = false, sortedby::Function = identity, presorted::Bool = false) where T
    sortedpoints = presorted ? points : sort(points; by = sortedby)
    pointslist = PairedLinkedList{T}(sortedpoints...)
    hull = PairedLinkedList{T}()
    addpartner!(hull, pointslist)
    h = MutableConvexHull{T}(hull, orientation, collinear, sortedby)
    monotonechain!(h)
    return h
end

