function insidehull(pointdata::T, h::MutableConvexHull{T}) where T
    length(h) == 1 && return pointdata == h.head.next.data

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
    length(h) == 1 && return pointdata == h.hull.head.next.data
    abovelower = false
    ccw = h.orientation == CCW
    for nextnode in h.hull.head.next.next
        prevnode = nextnode.prev
        if ccw ? (prevnode.data[1] <= pointdata[1] <= nextnode.data[1]) : (prevnode.data[1] >= pointdata[1] >= nextnode.data[1])
            yhull = linterp(pointdata[1], prevnode.data, nextnode.data)
            if pointdata[2] >= yhull
                abovelower = true
                break
            end
        end
    end
    return abovelower
end

function insidehull(pointdata::T, h::MutableUpperConvexHull{T}) where T
    length(h) == 1 && return pointdata == h.head.next.data
    belowupper = false
    ccw = h.orientation == CCW
    for nextnode in h.hull.head.next.next
        if ccw ? (prevnode.data[1] >= pointdata[1] >= nextnode.data[1]) : (prevnode.data[1] <= pointdata[1] <= nextnode.data[1])
            yhull = linterp(pointdata[1], prevnode.data, nextnode.data)
            if pointdata[2] <= yhull
                belowupper = true
                break
            end
        end
    end
    return abovelower && belowupper
end

insidehull(pointnode::PairedListNode, h::AbstractConvexHull) = insidehull(pointnode.data, h)