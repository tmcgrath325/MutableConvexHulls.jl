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

function jarvismarch!(pointslist::PairedLinkedList{T}, hull::Union{PairedLinkedList{T},Nothing}=nothing, stop::Union{PairedListNode{T},Nothing}=nothing, initedge=DOWN; orientation::HullOrientation=CCW, collinear::Bool=false, by::Function=identity) where T
    if isnothing(hull)
        # initialize the convex hull
        hull = PairedLinkedList{T}()
        addpartner!(pointslist, hull)
    end

    # handle the 0- and 1- point cases
    if length(pointslist) == 1 
        push!(hull, first(pointslist))
        addpartner!(tail(hull), head(pointslist))
    end
    length(pointslist) <= 1 && return hull

    if isempty(hull)
        # add first point on the hull
        f = node -> by(node.data)
        firstnode = orientation === CCW ? argmin(f, IteratingListNodes(pointslist)) : argmax(f, IteratingListNodes(pointslist))
        push!(hull, firstnode.data)
        addpartner!(head(hull), firstnode)
    end
    if isnothing(stop)
        stop = head(hull)
    end

    # use the appropriate check for determining a better option for the next point
    betterturn(args...) = collinear ? iscloserturn(!orientation, args...) : isfurtherturn(!orientation, args...)

    # perform jarvis march 
    counter = 0
    current = tail(hull).partner
    prevedge = initedge
    while counter == 0 || current !== stop.partner
        if counter > length(pointslist)
            throw(ErrorException("More points were added to the hull than exist in the provided list of points."))
        end
        counter += 1
        next = jarvissearch(current, prevedge, IteratingListNodes(pointslist), betterturn)
        if current == next
            throw(ErrorException("Jarvis March failed to progress."))
        end
        # stop when the first point on the hull has been reached (usually to return only the upper/lower convex hull)
        if next == stop
            break
        end
        # add the next node to the hull
        push!(hull, next.data)
        addpartner!(tail(hull), next)
        prevedge = next.data .- current.data
        current = next
    end

    return hull
end

function lower_jarvismarch!(pointslist::PairedLinkedList{T}; orientation::HullOrientation=CCW, collinear::Bool=false, by::Function=identity) where T
    # initialize the hull
    hull = PairedLinkedList{T}()
    addpartner!(pointslist, hull)

    # handle the 0- and 1-point cases
    length(pointslist) == 1 && push!(hull, first(pointslist))
    length(pointslist) <= 1 && return hull

    # select the appropriate starting and stopping nodes
    f = node -> by(node.data)
    firstnode = orientation === CCW ? argmin(f, IteratingListNodes(pointslist)) : argmax(f, IteratingListNodes(pointslist))
    stop = orientation === CW ? argmin(f, IteratingListNodes(pointslist)) : argmax(f, IteratingListNodes(pointslist))
    push!(hull, firstnode.data)
    addpartner!(head(hull), firstnode)

    # populate the hull via jarvis march
    jarvismarch!(pointslist, hull, stop, DOWN; orientation=orientation, collinear=collinear, by=by)

    # add the last node
    push!(hull, stop.data)
    addpartner!(tail(hull), stop)

    return hull
end

function upper_jarvismarch!(pointslist::PairedLinkedList{T}; orientation::HullOrientation=CCW, collinear::Bool=false, by::Function=identity) where T
    # initialize the hull
    hull = PairedLinkedList{T}()
    addpartner!(pointslist, hull)

    # handle the 0- and 1-point cases
    length(pointslist) == 1 && push!(hull, first(pointslist))
    length(pointslist) <= 1 && return hull

    # select the appropriate starting and stopping nodes
    f = node -> by(node.data)
    firstnode = orientation === CW ? argmin(f, IteratingListNodes(pointslist)) : argmax(f, IteratingListNodes(pointslist))
    stop = orientation === CCW ? argmin(f, IteratingListNodes(pointslist)) : argmax(f, IteratingListNodes(pointslist))
    push!(hull, firstnode.data)
    addpartner!(head(hull), firstnode)

    # populate the hull via jarvis march
    jarvismarch!(pointslist, hull, stop, UP; orientation=orientation, collinear=collinear, by=by)

    # add the last node
    push!(hull, stop.data)
    addpartner!(tail(hull), stop)

    return hull
end

function jarvismarch(points::AbstractVector{T}; orientation::HullOrientation=CCW, collinear::Bool=false, kwargs...) where T
    pointslist = PairedLinkedList{T}(points...)
    hull = jarvismarch!(pointslist; orientation=orientation, collinear=collinear, kwargs...)
    return MutableConvexHull{T}(hull, orientation, collinear)
end

function lower_jarvismarch(points::AbstractVector{T}; orientation::HullOrientation=CCW, collinear::Bool=false, kwargs...) where T
    pointslist = PairedLinkedList{T}(points...)
    hull = lower_jarvismarch!(pointslist; orientation=orientation, collinear=collinear, kwargs...)
    return MutableLowerConvexHull{T}(hull, orientation, collinear)
end

function upper_jarvismarch(points::AbstractVector{T}; orientation::HullOrientation=CCW, collinear::Bool=false, kwargs...) where T
    pointslist = PairedLinkedList{T}(points...)
    hull = upper_jarvismarch!(pointslist; orientation=orientation, collinear=collinear, kwargs...)
    return MutableUpperConvexHull{T}(hull, orientation, collinear)
end