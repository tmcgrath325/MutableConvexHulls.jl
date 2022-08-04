function jarvismarch!(pointslist::PairedLinkedList{T}, hull::Union{PairedLinkedList{T},Nothing}=nothing, stop::Union{PairedListNode{T},Nothing}=nothing; orientation::HullOrientation=CCW, colinear::Bool=false) where T
    if isnothing(hull)
        # initialize the convex hull
        hull = PairedLinkedList{T}()
        addpartner!(pointslist, hull)
    end

    # handle the 0- and 1-point cases
    length(pointslist) == 1 && push!(hull, first(pointslist))
    length(pointslist) <= 1 && return hull

    if isempty(hull)
        # add first point on the hull
        firstnode = getnode(pointslist, orientation === CCW ? argmin(pointslist) : argmax(pointslist))
        push!(hull, firstnode.data)
        addpartner!(hull.head.next, firstnode)
    end
    if isnothing(stop)
        stop = hull.head.next
    end

    # use the appropriate check for determining a better option for the next point
    betterturn(args...) = colinear ? closer_turn(!orientation, args...) : further_turn(!orientation, args...)

    # perform jarvis march 
    current = hull.tail.prev.partner
    firstpoint = pointslist.head.next
    while length(hull) < 1 || current !== first(hull)
        next = firstpoint === current ? firstpoint.next : firstpoint # avoid checking identical points
        for target in IteratingListNodes(pointslist)
            if target.data != current.data
                # update the next node if it presents a better turn
                if betterturn(current.data, next.data, target.data)
                    next = target
                end
            end
        end
        if current == next
            throw(ErrorException("Jarvis March failed to progress."))
        end
        # stop when the first point on the hull has been reached (usually to return only the upper/lower convex hull)
        if next == stop
            break
        end
        # add the next node to the hull
        push!(hull, next.data)
        addpartner!(hull.tail.prev, next)
        current = next
    end

    return hull
end

function lower_jarvismarch!(pointslist::PairedLinkedList{T}; orientation::HullOrientation=CCW, colinear::Bool=false) where T
    # initialize the hull
    hull = PairedLinkedList{T}()
    addpartner!(pointslist, hull)

    # handle the 0- and 1-point cases
    length(pointslist) == 1 && push!(hull, first(pointslist))
    length(pointslist) <= 1 && return hull

    # select the appropriate starting and stopping nodes
    firstnode = getnode(pointslist, orientation === CCW ? argmin(pointslist) : argmax(pointslist))
    stop = getnode(pointslist, orientation === CW ? argmin(pointslist) : argmax(pointslist))
    push!(hull, firstnode.data)
    addpartner!(hull.head.next, firstnode)

    # populate the hull via jarvis march
    jarvismarch!(pointslist, hull, stop; orientation=orientation, colinear=colinear)

    # add the last node
    push!(hull, stop.data)
    addpartner!(hull.tail.prev, stop)

    return hull
end

function upper_jarvismarch!(pointslist::PairedLinkedList{T}; orientation::HullOrientation=CCW, colinear::Bool=false) where T
    # initialize the hull
    hull = PairedLinkedList{T}()
    addpartner!(pointslist, hull)

    # handle the 0- and 1-point cases
    length(pointslist) == 1 && push!(hull, first(pointslist))
    length(pointslist) <= 1 && return hull

    # select the appropriate starting and stopping nodes
    firstnode = getnode(pointslist, orientation === CW ? argmin(pointslist) : argmax(pointslist))
    stop = getnode(pointslist, orientation === CCW ? argmin(pointslist) : argmax(pointslist))
    push!(hull, firstnode.data)
    addpartner!(hull.head.next, firstnode)

    # populate the hull via jarvis march
    jarvismarch!(pointslist, hull, stop; orientation=orientation, colinear=colinear)

    # add the last node
    push!(hull, stop.data)
    addpartner!(hull.tail.prev, stop)

    return hull
end

function jarvismarch(points::AbstractVector{T}; orientation::HullOrientation=CCW, colinear::Bool=false) where T
    pointslist = PairedLinkedList{T}(points...)
    hull = jarvismarch!(pointslist; orientation=orientation, colinear=colinear)
    return MutableConvexHull{T}(hull, orientation, colinear)
end

function lower_jarvismarch(points::AbstractVector{T}; orientation::HullOrientation=CCW, colinear::Bool=false) where T
    pointslist = PairedLinkedList{T}(points...)
    hull = lower_jarvismarch!(pointslist; orientation=orientation, colinear=colinear)
    return MutableLowerConvexHull{T}(hull, orientation, colinear)
end

function upper_jarvismarch(points::AbstractVector{T}; orientation::HullOrientation=CCW, colinear::Bool=false) where T
    pointslist = PairedLinkedList{T}(points...)
    hull = upper_jarvismarch!(pointslist; orientation=orientation, colinear=colinear)
    return MutableUpperConvexHull{T}(hull, orientation, colinear)
end