abstract type AbstractConvexHull{T} end

mutable struct MutableConvexHull{T, F<:Function} <: AbstractConvexHull{T}
    hull::HullList{T,F}
    points::PointList{T,HullList{T,F},HullNode{T,HullList{T,F},F},F}
    orientation::HullOrientation
    collinear::Bool
    sortedby::F
end
function MutableConvexHull{T,F}(orientation::HullOrientation=CCW, collinear::Bool=false, sortedby::F=identity) where {T,F}
    points = PointList{T}(;sortedby=sortedby)
    hull = HullList{T,F}()
    addtarget!(hull,points)
    return MutableConvexHull{T,F}(hull, points, orientation, collinear, sortedby)
end
MutableConvexHull{T}(orientation=CCW, collienar=false, sortedby::F=identity) where {T,F} = MutableConvexHull{T,F}(orientation,collienar,sortedby)

mutable struct MutableLowerConvexHull{T, F<:Function} <: AbstractConvexHull{T}
    hull::HullList{T,F}
    points::PointList{T,HullList{T,F},HullNode{T,HullList{T,F},F},F}
    orientation::HullOrientation
    collinear::Bool
    sortedby::F
end
function MutableLowerConvexHull{T,F}(orientation::HullOrientation=CCW, collinear::Bool=false, sortedby::F=identity) where {T,F}
    points = PointList{T}(;sortedby=sortedby)
    hull = HullList{T,F}()
    addtarget!(hull,points)
    return MutableLowerConvexHull{T,F}(hull, points, orientation, collinear, sortedby)
end
MutableLowerConvexHull{T}(orientation=CCW, collienar=false, sortedby::F=identity) where {T,F} = MutableLowerConvexHull{T,F}(orientation,collienar,sortedby)

mutable struct MutableUpperConvexHull{T, F<:Function} <: AbstractConvexHull{T}
    hull::HullList{T,F}
    points::PointList{T,HullList{T,F},HullNode{T,HullList{T,F},F},F}
    orientation::HullOrientation
    collinear::Bool
    sortedby::F
end
function MutableUpperConvexHull{T,F}(orientation::HullOrientation=CCW, collinear::Bool=false, sortedby::F=identity) where {T,F}
    points = PointList{T}(;sortedby=sortedby)
    hull = HullList{T,F}()
    addtarget!(hull,points)
    return MutableUpperConvexHull{T,F}(hull, points, orientation, collinear, sortedby)
end
MutableUpperConvexHull{T}(orientation=CCW, collienar=false, sortedby::F=identity) where {T,F} = MutableUpperConvexHull{T,F}(orientation,collienar, sortedby)

Base.isempty(h::AbstractConvexHull) = isempty(h.hull)
Base.length(h::AbstractConvexHull) = length(h.hull)
Base.eltype(h::AbstractConvexHull) = eltype(h.hull)

Base.:(==)(h1::AbstractConvexHull, h2::AbstractConvexHull) = h1.hull == h2.hull

# function Base.copy!(h2::H, h::H) where {T,H<:AbstractConvexHull{T}}
#     copy!(h2.hull, h.hull)
#     h2.points = h2.hull.target
#     h2.orientation = h.orientation
#     h2.collinear = h.collinear
#     h2.sortedby = h.sortedby
#     return h2
# end
# function Base.copy(h::H) where {T,H<:AbstractConvexHull{T}}
#     copiedhull = copy(h.hull)
#     return H(copiedhull, copiedhull.partner, h.orientation, h.collinear, h.sortedby, h.issorted)
# end

Base.empty!(h::AbstractConvexHull) = empty!(h.hull)
Base.empty(h::H) where H <: AbstractConvexHull = H(h.orientation, h.collinear, h.sortedby)

# Iterating a convex hull returns the data contained in the nodes of its hull list
Base.iterate(h::AbstractConvexHull) = iterate(h, h.hull.head.next)
Base.iterate(h::AbstractConvexHull, node::HullNode) = iterate(h.hull, node)

struct HullNodeIterator{S<:HullNode}
    start::S
    rev::Bool
end
"""
    HullNodeIterator(start [, rev])

Returns an iterator over the nodes in the linked list representing the convex hull, starting at the specified node `start`.

If `rev` is `true`, the iterator will advance toward the head of the list.
Otherwise, it will advance toward the tail of the list.
"""
function HullNodeIterator(start::S; rev::Bool = false) where {S<:HullNode}
    return HullNodeIterator{S}(start, rev)
end

"""
    HullNodeIterator(hull [, rev])

Returns an iterator over the nodes in the linked list representing the convex hull.

If `rev` is `true`, the iterator will start at the tail of the list and advance toward the head.
Otherwise, it will start at the head of the list and advance toward the tail.
"""
HullNodeIterator(h::AbstractConvexHull; rev::Bool = false) = HullNodeIterator{nodetype(h.hull)}(rev ? h.hull.tail.prev : h.hull.head.next, rev)

Base.iterate(iter::HullNodeIterator) = iterate(iter, iter.start)
Base.iterate(iter::HullNodeIterator{S}, node::S) where {S<:HullNode} = iter.rev ? (athead(node) ? nothing : (node, node.prev)) :
                                                                                  (attail(node) ? nothing : (node, node.next))
Base.IteratorSize(::HullNodeIterator) = Base.SizeUnknown()


struct PointNodeIterator{S<:PointNode}
    start::S
    rev::Bool
end
"""
    PointNodeIterator(start [, rev])

Returns an iterator over the nodes in the linked list representing the points contained by the convex hull,
starting at the specified node `start`.

If `rev` is `true`, the iterator will advance toward the head of the list.
Otherwise, it will advance toward the tail of the list.
"""
function PointNodeIterator(start::S; rev::Bool = false) where S
    return PointNodeIterator{S}(start, rev)
end

"""
    PointNodeIterator(hull [, rev])

Returns an iterator over the nodes in the linked list representing the points contained by the convex hull.

If `rev` is `true`, the iterator will start at the tail of the list and advance toward the head.
Otherwise, it will start at the head of the list and advance toward the tail.
"""
PointNodeIterator(h::AbstractConvexHull; rev::Bool = false) = PointNodeIterator{nodetype(h.points)}(rev ? h.points.tail.prev : h.points.head.next, rev)

Base.iterate(iter::PointNodeIterator) = iterate(iter, iter.start)
Base.iterate(iter::PointNodeIterator{S}, node::S) where S = iter.rev ? (athead(node) ? nothing : (node, node.prev)) :
                                                                       (attail(node) ? nothing : (node, node.next))
Base.IteratorSize(::PointNodeIterator) = Base.SizeUnknown()


# struct BracketedIteratorState{T}
#     node::PointNode{T}
#     finished::Bool
# end

# struct BracketedPointNodeIterator{T}
#     start::PointNode{T}
#     hullstart::HullNode{T}
#     hullend::HullNode{T}
#     rev::Bool
# end
# function BracketedPointNodeIterator(start::PointNode{T}, hullstart::HullNode{T}, hullend::HullNode{T}; rev::Bool = false) where T
#     hullstart.wrapped.list !== hullend.wrapped.list && throw(ArgumentError("The hull nodes must belong to the same convex hull"))
#     start.wrapped.list.target !== hullstart.wrapped.list && throw(ArgumentError("The starting point node must belong to the same convex hull as the hull points"))
#     return BracketedPointNodeIterator{T}(start, hullstart, hullend, rev)
# end

# Base.iterate(iter::BracketedPointNodeIterator{T}) where T = iterate(iter, BracketedIteratorState{T}(iter.start, iter.start.wrapped.list.len == 0))
# function Base.iterate(iter::BracketedPointNodeIterator{T}, state::BracketedIteratorState{T}) where T
#     state.finished && return nothing
#     if iter.rev
#         prevnode = athead(state.node.wrapped.prev) ? tail(state.node.wrapped.list) : state.node.wrapped.prev
#         return (state.node.wrapped, BracketedIteratorState(PointNode{T}(prevnode), state.node.wrapped.target === iter.hullstart.wrapped))
#     else
#         nextnode = attail(state.node.wrapped.next) ? head(state.node.wrapped.list) : state.node.wrapped.next
#         return (state.node.wrapped, BracketedIteratorState(PointNode{T}(nextnode), state.node.wrapped.target === iter.hullend.wrapped))
#     end
# end
# Base.IteratorSize(::BracketedPointNodeIterator) = Base.SizeUnknown()


function addpoint!(h::AbstractConvexHull{T}, point::T) where T
    # handle the case when the hull is initially empty
    if length(h) == 0
        push!(h.hull, point)
        push!(h.points, point)
        addtarget!(tail(h.hull), tail(h.points))
        return h
    end
    push!(h.points, point)
    if !insidehull(point, h)    # if the new point is outside the hull, update the convex hull
        monotonechain!(h)
    end
    return h
end

function mergepoints!(h::MutableConvexHull{T}, points::AbstractVector{T}) where T
    h2 = monotonechain(points; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
    mergehulls!(h,h2)
    return h
end
function mergepoints!(h::MutableLowerConvexHull{T}, points::AbstractVector{T}) where T
    h2 = lower_monotonechain(points; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
    mergehulls!(h,h2)
    return h
end
function mergepoints!(h::MutableUpperConvexHull{T}, points::AbstractVector{T}) where T
    h2 = upper_monotonechain(points; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
    mergehulls!(h,h2)
    return h
end

function removepoint!(h::AbstractConvexHull{T}, node::HullNode{T}) where T
    node.list !== h.hull && throw(ArgumentError("The specified node must belong to the provided convex hull"))
    start = node.prev.target
    stop = node.next.target
    deletenode!(node.target)
    deletenode!(node)
    start = start.list === h.points ? start : firstpoint(h)
    stop = stop.list === h.points ? stop : lastpoint(h)
    monotonechain!(h, start, stop)
    return h
end

function removepoint!(h::AbstractConvexHull{T}, node::PointNode{T}) where T
    node.list !== h.points && throw(ArgumentError("The specified node must belong to the provided convex hull"))
    if hastarget(node)
        start = node.target.prev.target
        stop = node.target.next.target
        deletenode!(node.target)
        deletenode!(node)
        start = start.list === h.points ? start : firstpoint(h)
        stop = stop.list === h.points ? stop : lastpoint(h)
        monotonechain!(h, start, stop)
    else
        deletenode!(node)
    end
    return h
end

function Base.show(io::IO, h::AbstractConvexHull)
    print(io, typeof(h), '(')
    join(io, h, ", ")
    print(io, ')')
end
