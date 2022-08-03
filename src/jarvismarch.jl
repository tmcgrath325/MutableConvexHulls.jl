Base.keys(l::PairedLinkedLists.AbstractLinkedList) = LinearIndices(1:l.len)

function jarvismarch!(pointslist::PairedLinkedList{T}, hull::Union{PairedLinkedList{T},Nothing}=nothing, stop::Union{PairedListNode{T},Nothing}=nothing; direction::Type{<:Union{CCW,CW}}=CCW, colinear::Bool=false) where T
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
        firstnode = getnode(pointslist, direction === CCW ? argmin(pointslist) : argmax(pointslist))
        push!(hull, firstnode.data)
        addpartner!(hull.head.next, firstnode)
    end
    if isnothing(stop)
        stop = hull.head.next
    end
    # use the appropriate check for determining a better option for the next point
    betterturn = direction === CCW ? 
        (colinear ? misaligned_right_turn : aligned_right_turn) :
        (colinear ? misaligned_left_turn : aligned_left_turn)

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
        # stop when the first point on the hull has been reached
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

function lower_jarvismarch!(pointslist::PairedLinkedList{T}; direction::Type{<:Union{CCW,CW}}=CCW, colinear::Bool=false) where T
    # initialize the hull
    hull = PairedLinkedList{T}()
    addpartner!(pointslist, hull)    
    # handle the 0- and 1-point cases
    length(pointslist) == 1 && push!(hull, first(pointslist))
    length(pointslist) <= 1 && return hull
    # select the appropriate starting and stopping nodes
    firstnode = getnode(pointslist, direction === CCW ? argmin(pointslist) : argmax(pointslist))
    stop = getnode(pointslist, direction === CW ? argmin(pointslist) : argmax(pointslist))
    push!(hull, firstnode.data)
    addpartner!(hull.head.next, firstnode)
    # populate the hull via jarvis march
    jarvismarch!(pointslist, hull, stop; direction=direction, colinear=colinear)
    # add the last node
    push!(hull, stop.data)
    addpartner!(hull.tail.prev, stop)
    return hull
end

function upper_jarvismarch!(pointslist::PairedLinkedList{T}; direction::Type{<:Union{CCW,CW}}=CCW, colinear::Bool=false) where T
    # initialize the hull
    hull = PairedLinkedList{T}()
    addpartner!(pointslist, hull)
    # handle the 0- and 1-point cases
    length(pointslist) == 1 && push!(hull, first(pointslist))
    length(pointslist) <= 1 && return hull
    # select the appropriate starting and stopping nodes
    firstnode = getnode(pointslist, direction === CW ? argmin(pointslist) : argmax(pointslist))
    stop = getnode(pointslist, direction === CCW ? argmin(pointslist) : argmax(pointslist))
    push!(hull, firstnode.data)
    addpartner!(hull.head.next, firstnode)
    # populate the hull via jarvis march
    jarvismarch!(pointslist, hull, stop; direction=direction, colinear=colinear)
    # add the last node
    push!(hull, stop.data)
    addpartner!(hull.tail.prev, stop)
    return hull
end
