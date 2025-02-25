function jarvissearch(query::AbstractNode{T}, prevedge, pointsnodes, betterturn::Function) where T
    # isempty(pointsnodes) && throw(ArgumentError("At least 1 point must be provided."))
    firstpoint, iter = Iterators.peel(pointsnodes)
    isempty(iter) && return firstpoint
    next = firstpoint === query ? first(iter) : firstpoint # avoid checking identical points
    for target in pointsnodes
        if !coordsareequal(target.data, query.data)
            # update the next node if it presents a better turn
            if betterturn(prevedge, query.data, next.data, target.data)
                next = target
            end
        end
    end
    return next
end

getiterator(pl) = typeof(pl) <: AbstractList ? ListNodeIterator(pl) : Iterators.Stateful(pl)

"""
    jarvismarch!(hull [, initedge, stop])

Determine the convex hull of the points contained in the provided `hull.points` using the [Jarvis march algorithm](https://en.wikipedia.org/wiki/Gift_wrapping_algorithm). 
Each node in the list should contain a two-dimensional point, and the nodes are assumed to be sorted 
(e.g. by lowest "x" value and by lowest "y" in case of ties, though some other sorting methods may produce valid results).
"""
function jarvismarch!(hull::Union{AbstractConvexHull{T},AbstractList{T}}, pointslist, collinear, orientation, initedge, stop::Union{PointNode{T},HullNode{T},Nothing}=nothing) where T
    # pointslist = h.hull.target   
    @assert length(hull) > 0
    if isnothing(stop)
        stop = head(hull)
    end
    stopdata = stop.data

    # use the appropriate check for determining a better option for the next point
    betterturn(prevedge,o,a,b) = collinear ? iscloserturn(!orientation,prevedge,o,a,b) : isfurtherturn(!orientation,prevedge,o,a,b)

    # perform jarvis march 
    counter = 0
    current = head(hull).target
    prevedge = initedge
    while counter == 0 || current !== stop.target
        if counter > length(pointslist)
            throw(ErrorException("More points were added to the hull than exist in the list of candidate points."))
        end
        counter += 1
        next = jarvissearch(current, prevedge, getiterator(pointslist), betterturn)
        if coordsareequal(current.data, next.data)
            if length(hull) == 1
                return hull
            else
                throw(ErrorException("Jarvis March failed to progress."))
            end
        end
        next == stop && break                                                       # the stopping point has been reached
        (typeof(pointslist) <: AbstractPointList{T}) && hastarget(next) && break    # the next node to be added already is part of the hull
        # add the next node to the hull
        push!(hull, next.data)
        addtarget!(tail(hull), next)
        prevedge = next.data .- current.data
        current = next
    end

    if length(hull) > 1
        if coordsareequal(head(hull).data, tail(hull).data) 
            deletenode!(tail(hull))
        end
    end

    return hull
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
    jarvismarch!(h.hull, h.hull.target, h.collinear, h.orientation, DOWN)

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

    # handle duplicate point cases
    coordsareequal(firstnode.data, stop.data) && return h

    # populate the hull via jarvis march
    jarvismarch!(h.hull, h.hull.target, h.collinear, h.orientation, DOWN, stop)

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

    # handle duplicate point cases
    coordsareequal(firstnode.data, stop.data) && return h

    # populate the hull via jarvis march
    jarvismarch!(h.hull, h.hull.target, h.collinear, h.orientation, UP, stop)

    # add the last node
    push!(h.hull, stop.data)
    addtarget!(tail(h.hull), stop)

    return h
end

"""
    h = jarvismarch(points [; orientation, collinear, sortedby])

Return the convex hull generated from the provided `points`.

`orientation` specifies whether the points along the convex hull are ordered clockwise `CW`, or counterclockwise `CCW`, and defaults to `CCW`.

`collinear` specifies whether collinear points are allowed along the surface of the convex hull, and defaults to `false`.

`sortedby` specifies a function to apply to points prior to sorting, and defaults to `identity` (resulting in default sorting behavior).
"""
function jarvismarch(points::AbstractVector{T}; orientation::HullOrientation=CCW, collinear::Bool=false, sortedby::Function=identity) where T
    pointslist = PointList{T}(;sortedby=sortedby)
    hull = HullList{T,typeof(sortedby)}()
    addtarget!(hull, pointslist)
    !isempty(points) && push!(pointslist, points...)
    h = MutableConvexHull{T, typeof(sortedby)}(hull, pointslist, orientation, collinear, sortedby)
    jarvismarch!(h)
    return h
end

"""
    lh = lower_jarvismarch(points [; orientation, collinear, sortedby])

Return the lower convex hull generated from the provided `points`.

`orientation` specifies whether the points along the convex hull are ordered clockwise `CW`, or counterclockwise `CCW`, and defaults to `CCW`.

`collinear` specifies whether collinear points are allowed along the surface of the convex hull, and defaults to `false`.

`sortedby` specifies a function to apply to points prior to sorting, and defaults to `identity` (resulting in default sorting behavior).
"""
function lower_jarvismarch(points::AbstractVector{T}; orientation::HullOrientation=CCW, collinear::Bool=false, sortedby::Function=identity) where T
    pointslist = PointList{T}(;sortedby=sortedby)
    hull = HullList{T,typeof(sortedby)}()
    addtarget!(hull, pointslist)
    !isempty(points) && push!(pointslist, points...)
    h = MutableLowerConvexHull{T, typeof(sortedby)}(hull, pointslist, orientation, collinear, sortedby)
    jarvismarch!(h)
    return h
end

"""
    uh = upper_jarvismarch(points [; orientation, collinear, sortedby])

Return the upper convex hull generated from the provided `points`.

`orientation` specifies whether the points along the convex hull are ordered clockwise `CW`, or counterclockwise `CCW`, and defaults to `CCW`.

`collinear` specifies whether collinear points are allowed along the surface of the convex hull, and defaults to `false`.

`sortedby` specifies a function to apply to points prior to sorting, and defaults to `identity` (resulting in default sorting behavior).
"""
function upper_jarvismarch(points::AbstractVector{T}; orientation::HullOrientation=CCW, collinear::Bool=false, sortedby::Function=identity) where T
    pointslist = PointList{T}(;sortedby=sortedby)
    hull = HullList{T,typeof(sortedby)}()
    addtarget!(hull, pointslist)
    !isempty(points) && push!(pointslist, points...)
    h = MutableUpperConvexHull{T, typeof(sortedby)}(hull, pointslist, orientation, collinear, sortedby)
    jarvismarch!(h)
    return h
end