"""
    insidehull(data, hull) -> Bool

Return `true` if the data lies within the interior of `hull` and `false` otherwise.
"""
function insidehull(pointdata::T, h::MutableConvexHull{T}) where T
    length(h) == 0 && return false
    length(h) == 1 && return pointdata == h.hull.head.next.data

    abovelower = false
    belowupper = false
    alreadycheckedlower = false
    ccw = h.orientation == CCW

    for prevnode in ListNodeIterator(h.hull)
        nextnode = attail(prevnode.next) ? head(h.hull) : prevnode.next
        nextnextnode = attail(nextnode.next) ? head(h.hull) : nextnode.next

        # the lower hull should always be listed first
        if !alreadycheckedlower
            # check if the lower hull has been passed through
            if ccw ? (prevnode.data[1] > nextnode.data[1]) : (prevnode.data[1] < nextnode.data[1])
                alreadycheckedlower = true
            elseif !abovelower
                # if the point lies along the extreme left or right edge of the entire hull...
                if pointdata[1] == prevnode.data[1] == nextnode.data[1] 
                    if h.collinear 
                        abovelower = (pointdata == nextnode.data || pointdata == prevnode.data) 
                    else
                        abovelower = nextnode.data[2] >= prevnode.data[2] ? prevnode.data[2] <= pointdata[2] <= nextnode.data[2] :
                                                                            prevnode.data[2] >= pointdata[2] >= nextnode.data[2]
                    end
                    abovelower && return true
                    nextnode.data != nextnextnode.data && return false
                # if the point is even with the previous node...
                elseif pointdata[1] == prevnode.data[1]
                    abovelower = h.collinear ? pointdata[2] > prevnode.data[2] : pointdata[2] >= prevnode.data[2]
                # if the point is even with the next node...
                elseif pointdata[1] == nextnode.data[1] !== nextnextnode.data[1]
                    abovelower = h.collinear ? pointdata[2] > nextnode.data[2] : pointdata[2] >= nextnode.data[2]
                # if the point is in between the previous and next nodes...
                elseif ccw ? (prevnode.data[1] < pointdata[1] < nextnode.data[1]) : (prevnode.data[1] > pointdata[1] > nextnode.data[1])
                    yhull = linterp(pointdata[1], prevnode.data, nextnode.data)
                    if h.collinear ? pointdata[2] > yhull : pointdata[2] >= yhull
                        abovelower = true
                    end
                end
            end
        end
        if alreadycheckedlower
            !abovelower && return false # we can stop if we already know the point is outside the hull
            # if the point lies along the extreme left or right edge of the entire hull...
            if pointdata[1] == prevnode.data[1] == nextnode.data[1]
                if h.collinear 
                    belowupper = (pointdata == nextnode.data || pointdata == prevnode.data) 
                else
                    belowupper = nextnode.data[2] >= prevnode.data[2] ? prevnode.data[2] <= pointdata[2] <= nextnode.data[2] :
                                                                        prevnode.data[2] >= pointdata[2] >= nextnode.data[2]
                end
                belowupper && return true # we can stop if we know the hull is within both the upper and lower hulls
            # if the point is even with the previous node...
            elseif pointdata[1] == prevnode.data[1]
                belowupper = h.collinear ? pointdata[2] < prevnode.data[2] : pointdata[2] <= prevnode.data[2]
            # if the point is even with the next node...
            elseif pointdata[1] == nextnode.data[1]
                belowupper = h.collinear ? pointdata[2] < nextnode.data[2] : pointdata[2] <= nextnode.data[2]
            # if the point is in between the previous and next nodes...
            elseif ccw ? (prevnode.data[1] >= pointdata[1] >= nextnode.data[1]) : (prevnode.data[1] <= pointdata[1] <= nextnode.data[1])
                yhull = linterp(pointdata[1], prevnode.data, nextnode.data)
                if h.collinear ? pointdata[2] < yhull : pointdata[2] <= yhull
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

    # handle cases where the new point is tied with an extreme point with respect to the first coordinate
    hhead, htail = head(h.hull), tail(h.hull)
    if pointdata[1] == hhead.data[1]
        lastleftnode = hhead
        for leftnode in ListNodeIterator(h.hull)
            leftnode.data == pointdata && return true
            leftnode.data[1] != hhead.data[1] && break
            lastleftnode = leftnode
        end
        if ccw 
            return pointdata[2] >= lastleftnode.data[2] && (!h.collinear || lastleftnode === hhead) && h.sortedby(pointdata) >= h.sortedby(hhead.data)
        else
            return pointdata[2] >= lastleftnode.data[2] && (!h.collinear || lastleftnode === hhead) && h.sortedby(pointdata) <= h.sortedby(hhead.data)
        end
    end
    if pointdata[1] == htail.data[1]
        firstrightnode = htail
        for rightnode in ListNodeIterator(h.hull; rev=true)
            rightnode.data == pointdata && return true
            rightnode.data[1] != htail.data[1] && break
            firstrightnode = rightnode
        end
        if ccw 
            return pointdata[2] >= firstrightnode.data[2] && (!h.collinear || firstrightnode === htail) && h.sortedby(pointdata) <= h.sortedby(htail.data)
        else
            return pointdata[2] >= firstrightnode.data[2] && (!h.collinear || firstrightnode === htail) && h.sortedby(pointdata) >= h.sortedby(htail.data)
        end
    end

    # all other cases
    for nextnode in h.hull.head.next.next
        prevnode = nextnode.prev
        if ccw ? (prevnode.data[1] <= pointdata[1] <= nextnode.data[1]) : (prevnode.data[1] >= pointdata[1] >= nextnode.data[1])
            yhull = linterp(pointdata[1], prevnode.data, nextnode.data)
            if h.collinear ? pointdata[2] > yhull : pointdata[2] >= yhull
                return true
            end
        end
    end
    return false
end

function insidehull(pointdata::T, h::MutableUpperConvexHull{T}) where T
    length(h) == 0 && return false
    length(h) == 1 && return pointdata == h.hull.head.next.data
    ccw = h.orientation == CCW

    # handle cases where the new point is tied with an extreme point with respect to the first coordinate
    hhead, htail = head(h.hull), tail(h.hull)
    if pointdata[1] == hhead.data[1]
        lastleftnode = hhead
        for leftnode in ListNodeIterator(h.hull)
            leftnode.data == pointdata && return true
            leftnode.data[1] != hhead.data[1] && break
            lastleftnode = leftnode
        end
        if ccw 
            return pointdata[2] <= lastleftnode.data[2] && (!h.collinear || lastleftnode === hhead) && h.sortedby(pointdata) <= h.sortedby(hhead.data)
        else
            return pointdata[2] <= lastleftnode.data[2] && (!h.collinear || lastleftnode === hhead) && h.sortedby(pointdata) >= h.sortedby(hhead.data)
        end
    end
    if pointdata[1] == htail.data[1]
        firstrightnode = htail
        for rightnode in ListNodeIterator(h.hull; rev=true)
            rightnode.data == pointdata && return true
            rightnode.data[1] != htail.data[1] && break
            firstrightnode = rightnode
        end
        if ccw 
            return pointdata[2] <= firstrightnode.data[2] && (!h.collinear || firstrightnode === htail) && h.sortedby(pointdata) >= h.sortedby(htail.data)
        else
            return pointdata[2] <= firstrightnode.data[2] && (!h.collinear || firstrightnode === htail) && h.sortedby(pointdata) <= h.sortedby(htail.data)
        end
    end

    # all other cases
    for nextnode in h.hull.head.next.next
        prevnode = nextnode.prev
        if ccw ? (prevnode.data[1] >= pointdata[1] >= nextnode.data[1]) : (prevnode.data[1] <= pointdata[1] <= nextnode.data[1])
            yhull = linterp(pointdata[1], prevnode.data, nextnode.data)
            if h.collinear ? pointdata[2] < yhull : pointdata[2] <= yhull
                return true
            end
        end
    end
    return false
end

insidehull(pointnode::PairedListNode, h::AbstractConvexHull) = insidehull(pointnode.data, h)