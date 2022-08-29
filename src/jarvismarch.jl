function jarvissearch(query::AbstractNode{T}, prevedge, pointsnodes, betterturn::Function) where T
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

function jarvismarch!(h::AbstractConvexHull{T}, initedge, stop::Union{PointNode{T},Nothing}=nothing) where T
    pointslist = h.hull.target
    hull = h.hull
    if isnothing(stop)
        stop = head(hull)
    end

    # use the appropriate check for determining a better option for the next point
    betterturn(args...) = h.collinear ? iscloserturn(!h.orientation, args...) : isfurtherturn(!h.orientation, args...)

    # perform jarvis march 
    counter = 0
    current = tail(hull).target
    prevedge = initedge
    while counter == 0 || current !== stop.target
        if counter > length(pointslist)
            throw(ErrorException("More points were added to the hull than exist in the list of candidate points."))
        end
        counter += 1
        next = jarvissearch(current, prevedge, ListNodeIterator(pointslist), betterturn)
        if current == next
            throw(ErrorException("Jarvis March failed to progress."))
        end
        # stop early when ...  
        next == stop && break           # the stopping point has been reached
        hastarget(next) && break       # the next node to be added already is part of the hull
        # add the next node to the hull
        push!(hull, next.data)
        addtarget!(tail(hull), next)
        prevedge = next.data .- current.data
        current = next
    end

    return h
end

function jarvismarch!(h::MutableConvexHull)
    # empty the hull and
    empty!(h.hull) # TODO: can we avoid starting over?

    # handle the 0- and 1-point cases
    if length(h.hull.target) == 1 
        push!(h.hull, first(h.hull.target))
        addtarget!(head(h.hull), head(h.hull.target))
    end
    length(h.hull.target) <= 1 && return h

    # select the appropriate starting node
    f = node -> h.sortedby(node.data)
    firstnode = h.orientation === CCW ? argmin(f, ListNodeIterator(h.hull.target)) : argmax(f, ListNodeIterator(h.hull.target))

    # reinitialize with the starting point
    push!(h.hull, firstnode.data)
    addtarget!(head(h.hull), firstnode)

    # populate the hull via jarvis march
    jarvismarch!(h, DOWN)

    return h
end

function jarvismarch!(h::MutableLowerConvexHull)
    # empty the hull and
    empty!(h.hull) # TODO: can we avoid starting over?

    # handle the 0- and 1-point cases
    if length(h.hull.target) == 1 
        push!(h.hull, first(h.hull.target))
        addtarget!(head(h.hull), head(h.hull.target))
    end
    length(h.hull.target) <= 1 && return h

    # select the appropriate starting and stopping nodes
    f = node -> h.sortedby(node.data)
    firstnode = h.orientation === CCW ? argmin(f, ListNodeIterator(h.hull.target)) : argmax(f, ListNodeIterator(h.hull.target))
    stop = h.orientation === CW ? argmin(f, ListNodeIterator(h.hull.target)) : argmax(f, ListNodeIterator(h.hull.target))

    # reinitialize with the starting point
    push!(h.hull, firstnode.data)
    addtarget!(head(h.hull), firstnode)

    # populate the hull via jarvis march
    jarvismarch!(h, DOWN, stop)

    # add the last node
    push!(h.hull, stop.data)
    addtarget!(tail(h.hull), stop)

    return h
end

function jarvismarch!(h::MutableUpperConvexHull)
    # empty the hull and
    empty!(h.hull) # TODO: can we avoid starting over?

    # handle the 0- and 1-point cases
    if length(h.hull.target) == 1 
        push!(h.hull, first(h.hull.target))
        addtarget!(head(h.hull), head(h.hull.target))
    end
    length(h.hull.target) <= 1 && return h.hull

    # select the appropriate starting and stopping nodes
    f = node -> h.sortedby(node.data)
    firstnode = h.orientation === CW ? argmin(f, ListNodeIterator(h.hull.target)) : argmax(f, ListNodeIterator(h.hull.target))
    stop = h.orientation === CCW ? argmin(f, ListNodeIterator(h.hull.target)) : argmax(f, ListNodeIterator(h.hull.target))

    # reinitialize with the starting point
    push!(h.hull, firstnode.data)
    addtarget!(head(h.hull), firstnode)

    # populate the hull via jarvis march
    jarvismarch!(h, UP, stop)

    # add the last node
    push!(h.hull, stop.data)
    addtarget!(tail(h.hull), stop)

    return h
end

function jarvismarch(points::AbstractVector{T}; orientation::HullOrientation=CCW, collinear::Bool=false, sortedby::Function=identity) where T
    pointslist = PointList{T}(;sortedby=sortedby)
    hull = HullList{T,typeof(sortedby)}()
    addtarget!(hull, pointslist)
    !isempty(points) && push!(pointslist, points...)
    h = MutableConvexHull{T, typeof(sortedby)}(hull, pointslist, orientation, collinear, sortedby)
    jarvismarch!(h)
    return h
end

function lower_jarvismarch(points::AbstractVector{T}; orientation::HullOrientation=CCW, collinear::Bool=false, sortedby::Function=identity) where T
    pointslist = PointList{T}(;sortedby=sortedby)
    hull = HullList{T,typeof(sortedby)}()
    addtarget!(hull, pointslist)
    !isempty(points) && push!(pointslist, points...)
    h = MutableLowerConvexHull{T, typeof(sortedby)}(hull, pointslist, orientation, collinear, sortedby)
    jarvismarch!(h)
    return h
end

function upper_jarvismarch(points::AbstractVector{T}; orientation::HullOrientation=CCW, collinear::Bool=false, sortedby::Function=identity) where T
    pointslist = PointList{T}(;sortedby=sortedby)
    hull = HullList{T,typeof(sortedby)}()
    addtarget!(hull, pointslist)
    !isempty(points) && push!(pointslist, points...)
    h = MutableUpperConvexHull{T, typeof(sortedby)}(hull, pointslist, orientation, collinear, sortedby)
    jarvismarch!(h)
    return h
end