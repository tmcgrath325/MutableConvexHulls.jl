function jarvissearch(query::AbstractListNode{T}, prevedge, pointsnodes, betterturn::Function) where T
    firstpoint = first(pointsnodes)
    next = firstpoint === query ? firstpoint.next : firstpoint # avoid checking identical points
    for target in pointsnodes
        if target.data != query.data
            # update the next node if it presents a better turn
            if betterturn(prevedge, query.data, next.data, target.data)
                next = target
            end
        end
    end
    return next
end

function jarvismarch!(h::AbstractConvexHull, initedge, stop::Union{PairedListNode{T},Nothing}=nothing) where T
    pointslist = h.hull.partner
    hull = h.hull
    if isnothing(stop)
        stop = head(hull)
    end

    # use the appropriate check for determining a better option for the next point
    betterturn(args...) = h.collinear ? iscloserturn(!h.orientation, args...) : isfurtherturn(!h.orientation, args...)

    # perform jarvis march 
    counter = 0
    current = tail(hull).partner
    prevedge = initedge
    while counter == 0 || current !== stop.partner
        if counter > length(pointslist)
            throw(ErrorException("More points were added to the hull than exist in the provided list of points."))
        end
        counter += 1
        next = jarvissearch(current, prevedge, ListNodeIterator(pointslist), betterturn)
        if current == next
            throw(ErrorException("Jarvis March failed to progress."))
        end
        # stop early when ...  
        next == stop && break           # the stopping point has been reached
        haspartner(next) && break       # the next node to be added already is part of the hull
        # add the next node to the hull
        push!(hull, next.data)
        addpartner!(tail(hull), next)
        prevedge = next.data .- current.data
        current = next
    end

    return h
end

function jarvismarch!(h::MutableConvexHull) where T
    # empty the hull and
    empty!(h.hull) # TODO: can we avoid starting over?

    # handle the 0- and 1-point cases
    length(h.hull.partner) == 1 && push!(hull, first(h.hull.partner))
    length(h.hull.partner) <= 1 && return hull

    # select the appropriate starting node
    f = isnothing(h.sortedby) ? node -> node.data : node -> h.sortedby(node.data)
    firstnode = h.orientation === CCW ? argmin(f, ListNodeIterator(h.hull.partner)) : argmax(f, ListNodeIterator(h.hull.partner))

    # reinitialize with the starting point
    push!(h.hull, firstnode.data)
    addpartner!(head(h.hull), firstnode)

    # populate the hull via jarvis march
    jarvismarch!(h, DOWN)

    return h
end

function jarvismarch!(h::MutableLowerConvexHull) where T
    # empty the hull and
    empty!(h.hull) # TODO: can we avoid starting over?

    # handle the 0- and 1-point cases
    length(h.hull.partner) == 1 && push!(hull, first(h.hull.partner))
    length(h.hull.partner) <= 1 && return hull

    # select the appropriate starting and stopping nodes
    f = isnothing(h.sortedby) ? node -> node.data : node -> h.sortedby(node.data)
    firstnode = h.orientation === CCW ? argmin(f, ListNodeIterator(h.hull.partner)) : argmax(f, ListNodeIterator(h.hull.partner))
    stop = h.orientation === CW ? argmin(f, ListNodeIterator(h.hull.partner)) : argmax(f, ListNodeIterator(h.hull.partner))

    # reinitialize with the starting point
    push!(h.hull, firstnode.data)
    addpartner!(head(h.hull), firstnode)

    # populate the hull via jarvis march
    jarvismarch!(h, DOWN, stop)

    # add the last node
    push!(h.hull, stop.data)
    addpartner!(tail(h.hull), stop)

    return h
end

function jarvismarch!(h::MutableUpperConvexHull) where T
    # empty the hull and
    empty!(h.hull) # TODO: can we avoid starting over?

    # handle the 0- and 1-point cases
    length(h.hull.partner) == 1 && push!(hull, first(h.hull.partner))
    length(h.hull.partner) <= 1 && return hull

    # select the appropriate starting and stopping nodes
    f = isnothing(h.sortedby) ? node -> node.data : node -> h.sortedby(node.data)
    firstnode = h.orientation === CW ? argmin(f, ListNodeIterator(h.hull.partner)) : argmax(f, ListNodeIterator(h.hull.partner))
    stop = h.orientation === CCW ? argmin(f, ListNodeIterator(h.hull.partner)) : argmax(f, ListNodeIterator(h.hull.partner))

    # reinitialize with the starting point
    push!(h.hull, firstnode.data)
    addpartner!(head(h.hull), firstnode)

    # populate the hull via jarvis march
    jarvismarch!(h, UP, stop)

    # add the last node
    push!(h.hull, stop.data)
    addpartner!(tail(h.hull), stop)

    return h
end

function jarvismarch(points::AbstractVector{T}; orientation::HullOrientation=CCW, collinear::Bool=false, sortedby::Union{Function,Nothing}=nothing) where T
    pointslist = PairedLinkedList{T}(points...)
    hull = PairedLinkedList{T}()
    addpartner!(hull, pointslist)
    h = MutableConvexHull{T}(hull, orientation, collinear, sortedby)
    jarvismarch!(h)
    return h
end

function lower_jarvismarch(points::AbstractVector{T}; orientation::HullOrientation=CCW, collinear::Bool=false, sortedby::Union{Function,Nothing}=nothing) where T
    pointslist = PairedLinkedList{T}(points...)
    hull = PairedLinkedList{T}()
    addpartner!(hull, pointslist)
    h = MutableLowerConvexHull{T}(hull, orientation, collinear, sortedby)
    jarvismarch!(h)
    return h
end

function upper_jarvismarch(points::AbstractVector{T}; orientation::HullOrientation=CCW, collinear::Bool=false, sortedby::Union{Function,Nothing}=nothing) where T
    pointslist = PairedLinkedList{T}(points...)
    hull = PairedLinkedList{T}()
    addpartner!(hull, pointslist)
    h = MutableUpperConvexHull{T}(hull, orientation, collinear, sortedby)
    jarvismarch!(h)
    return h
end