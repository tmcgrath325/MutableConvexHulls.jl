function insidehull(pointdata::T, h::MutableConvexHull{T}) where T
    length(h) == 0 && return false
    length(h) == 1 && return pointdata == h.hull.head.next.data

    abovelower = false
    belowupper = false
    alreadycheckedlower = false
    ccw = h.orientation == CCW

    for nextnode in h.hull.head.next.next
        prevnode = nextnode.prev
        # the lower hull should always be listed first
        if ccw ? (prevnode.data[1] > nextnode.data[1]) : (prevnode.data[1] < nextnode.data[1])
            alreadycheckedlower = true
        end
        if !alreadycheckedlower
            if ccw ? (prevnode.data[1] <= pointdata[1] <= nextnode.data[1]) : (prevnode.data[1] >= pointdata[1] >= nextnode.data[1])
                yhull = linterp(pointdata[1], prevnode.data, nextnode.data)
                if pointdata[2] >= yhull
                    abovelower = true
                    alreadycheckedlower = true
                end
            end
        end
        if alreadycheckedlower
            !abovelower && return false
            if ccw ? (prevnode.data[1] >= pointdata[1] >= nextnode.data[1]) : (prevnode.data[1] <= pointdata[1] <= nextnode.data[1])
                yhull = linterp(pointdata[1], prevnode.data, nextnode.data)
                if pointdata[2] <= yhull
                    belowupper = true
                    break
                end
            end
        end
    end
    return abovelower && belowupper
end

function insidehull(pointdata::T, h::MutableLowerConvexHull{T}) where T
    length(h) == 0 && return false
    length(h) == 1 && return pointdata == h.hull.head.next.data
    ccw = h.orientation === CCW

    # handle cases where the new point is tied with an extreme point
    hhead, htail = head(h.hull), tail(h.hull)
    if pointdata[1] == hhead.data[1]
        lastleftnode = hhead
        for leftnode in ListNodeIterator(h.hull)
            leftnode.data[1] != hhead.data[1] && break
            lastleftnode = leftnode
        end
        if lastleftnode === hhead
            return ccw ? h.sortedby(hhead.data) <= h.sortedby(pointdata) : h.sortedby(hhead.data) >= h.sortedby(pointdata)
        else
            @show lastleftnode.data, hhead.data
            return lastleftnode.data[2] > hhead.data[2] ? (hhead.data[2] <= pointdata[2] <= lastleftnode.data[2]) : 
                                                          (hhead.data[2] >= pointdata[2] >= lastleftnode.data[2])
        end
    end
    if pointdata[1] == htail.data[1]
        firstrightnode = htail
        for rightnode in ListNodeIterator(h.hull; rev=true)
            rightnode.data[1] != htail.data[1] && break
            firstrightnode = rightnode
        end
        if firstrightnode === htail
            return ccw ? h.sortedby(pointdata) <= h.sortedby(htail.data) : h.sortedby(pointdata) >= h.sortedby(hhead.data)
        else
            return firstrightnode.data[2] < htail.data[2] ? (firstrightnode.data[2] <= pointdata[2] <= htail.data[2]) : 
                                                            (firstrightnode.data[2] >= pointdata[2] >= htail.data[2])
        end
    end

    # all other cases
    for nextnode in h.hull.head.next.next
        prevnode = nextnode.prev
        if ccw ? (prevnode.data[1] <= pointdata[1] <= nextnode.data[1]) : (prevnode.data[1] >= pointdata[1] >= nextnode.data[1])
            yhull = linterp(pointdata[1], prevnode.data, nextnode.data)
            if pointdata[2] >= yhull
                return true
            end
        end
    end
    return false
end

function insidehull(pointdata::T, h::MutableUpperConvexHull{T}) where T
    length(h) == 0 && return false
    length(h) == 1 && return pointdata == h.head.next.data
    ccw = h.orientation == CCW

    # handle cases where the new point is tied with an extreme point
    hhead, htail = head(h.hull), tail(h.hull)
    if pointdata[1] == hhead.data[1]
        lastleftnode = hhead
        for leftnode in ListNodeIterator(h.hull)
            leftnode.data[1] != hhead.data[1] && break
            lastleftnode = leftnode
        end
        @show lastleftnode.data
        if lastleftnode === hhead
            return ccw ? h.sortedby(hhead.data) >= h.sortedby(pointdata) : h.sortedby(hhead.data) <= h.sortedby(pointdata)
        else
            return lastleftnode.data[2] > hhead.data[2] ? (hhead.data[2] <= pointdata[2] <= lastleftnode.data[2]) : 
                                                          (hhead.data[2] >= pointdata[2] >= lastleftnode.data[2])
        end
    end
    if pointdata[1] == htail.data[1]
        firstrightnode = htail
        for rightnode in ListNodeIterator(h.hull; rev=true)
            rightnode.data[1] != htail.data[1] && break
            firstrightnode = rightnode
        end
        @show firstrightnode.data, htail.data
        if firstrightnode === htail
            @show "here"
            return ccw ? h.sortedby(pointdata) >= h.sortedby(htail.data) : h.sortedby(pointdata) <= h.sortedby(hhead.data)
        else
            return firstrightnode.data[2] < htail.data[2] ? (firstrightnode.data[2] <= pointdata[2] <= htail.data[2]) : 
                                                            (firstrightnode.data[2] >= pointdata[2] >= htail.data[2])
        end
    end

    # all other cases
    for nextnode in h.hull.head.next.next
        if ccw ? (prevnode.data[1] >= pointdata[1] >= nextnode.data[1]) : (prevnode.data[1] <= pointdata[1] <= nextnode.data[1])
            yhull = linterp(pointdata[1], prevnode.data, nextnode.data)
            if pointdata[2] <= yhull
                return true
            end
        end
    end
    return false
end

insidehull(pointnode::PairedListNode, h::AbstractConvexHull) = insidehull(pointnode.data, h)