abstract type AbstractConvexHull{T} end

mutable struct MutableConvexHull{T, F<:Function} <: AbstractConvexHull{T}
    hull::PairedLinkedList{T}
    orientation::HullOrientation
    collinear::Bool
    sortedby::F
    issorted::Bool
end
function MutableConvexHull{T}(orientation::HullOrientation=CCW, collinear::Bool=false, sortedby::F=identity, issorted::Bool=false) where {T,F}
    pointslist = PairedLinkedList{T}()
    hull = PairedLinkedList{T}()
    addpartner!(hull,pointslist)
    return MutableConvexHull{T,F}(hull, orientation, collinear, sortedby, issorted)
end
MutableConvexHull{T,F}(orientation, collienar, sortedby, issorted) where {T,F} = MutableConvexHull{T}(orientation,collienar, sortedby, issorted)

mutable struct MutableLowerConvexHull{T, F<:Function} <: AbstractConvexHull{T}
    hull::PairedLinkedList{T}
    orientation::HullOrientation
    collinear::Bool
    sortedby::F
    issorted::Bool
end
function MutableLowerConvexHull{T}(orientation::HullOrientation=CCW, collinear::Bool=false, sortedby::F=identity, issorted::Bool=false) where {T,F}
    pointslist = PairedLinkedList{T}()
    hull = PairedLinkedList{T}()
    addpartner!(hull,pointslist)
    return MutableLowerConvexHull{T,F}(hull, orientation, collinear, sortedby, issorted)
end
MutableLowerConvexHull{T,F}(orientation, collienar, sortedby, issorted) where {T,F} = MutableLowerConvexHull{T}(orientation,collienar, sortedby, issorted)

mutable struct MutableUpperConvexHull{T, F<:Function} <: AbstractConvexHull{T}
    hull::PairedLinkedList{T}
    orientation::HullOrientation
    collinear::Bool
    sortedby::F
    issorted::Bool
end
function MutableUpperConvexHull{T}(orientation::HullOrientation=CCW, collinear::Bool=false, sortedby::F=identity, issorted::Bool=false) where {T,F}
    pointslist = PairedLinkedList{T}()
    hull = PairedLinkedList{T}()
    addpartner!(hull,pointslist)
    return MutableUpperConvexHull{T,F}(hull, orientation, collinear, sortedby, issorted)
end
MutableUpperConvexHull{T,F}(orientation, collienar, sortedby, issorted) where {T,F} = MutableUpperConvexHull{T}(orientation,collienar, sortedby, issorted)

Base.isempty(h::AbstractConvexHull) = isempty(h.hull)
Base.length(h::AbstractConvexHull) = length(h.hull)
Base.eltype(h::AbstractConvexHull) = eltype(h.hull)

Base.:(==)(h1::AbstractConvexHull, h2::AbstractConvexHull) = h1.hull == h2.hull

function Base.copy!(h2::H, h::H) where {T,H<:AbstractConvexHull{T}}
    copy!(h2.hull, h.hull)
    h2.orientation = h.orientation
    h2.collinear = h.collinear
    h2.sortedby = h.sortedby
    h2.issorted = h.issorted
    return h2
end
function Base.copy(h::H) where {T,H<:AbstractConvexHull{T}}
    copiedhull = copy(h.hull)
    return H(copiedhull, h.orientation, h.collinear, h.sortedby, h.issorted)
end

# function Base.reverse!(h::AbstractConvexHull)
#     h.hull = reverse!(hull)
#     h.orientation = !h.orientation
#     h.sortedby = x -> -1 .* h.sortedby(x)
#     return h
# end
# Base.reverse(h::AbstractConvexHull) = reverse!(copy(h))

Base.empty!(h::AbstractConvexHull) = empty!(h.hull)
Base.empty(h::H) where H <: AbstractConvexHull = H(h.orientation, h.collinear, h.sortedby, h.issorted)

# Iterating a convex hull returns the data contained in the nodes of its hull list
Base.iterate(h::AbstractConvexHull) = iterate(h, h.hull.head.next)
Base.iterate(h::AbstractConvexHull, node::PairedListNode) = iterate(h.hull, node)

# these exist purely to attempt to prevent accidentally misuse of HullNodeIterator, PointNodeIterator, and BracketedPointNodeIterator
struct PointNode{T}
    wrapped::PairedListNode{T}
end
struct HullNode{T}
    wrapped::PairedListNode{T}
end


struct HullNodeIterator{T}
    start::HullNode{T}
    rev::Bool
end
"""
    HullNodeIterator(start [, rev])

Returns an iterator over the nodes in the linked list representing the convex hull, starting at the specified node `start`.

If `rev` is `true`, the iterator will advance toward the head of the list.
Otherwise, it will advance toward the tail of the list.
"""
function HullNodeIterator(start::S; rev::Bool = false) where {T, S<:HullNode{T}}
    return HullNodeIterator{T}(start, rev)
end

"""
    HullNodeIterator(hull [, rev])

Returns an iterator over the nodes in the linked list representing the convex hull.

If `rev` is `true`, the iterator will start at the tail of the list and advance toward the head.
Otherwise, it will start at the head of the list and advance toward the tail.
"""
HullNodeIterator(h::AbstractConvexHull{T}; rev::Bool = false) where T = HullNodeIterator(rev ? HullNode{T}(h.hull.tail.prev) : HullNode{T}(h.hull.head.next); rev = rev)

Base.iterate(iter::HullNodeIterator) = iterate(iter, iter.start)
Base.iterate(iter::HullNodeIterator{T}, node::S) where {T, S<:HullNode{T}} = iter.rev ? (athead(node.wrapped) ? nothing : (node.wrapped, HullNode{T}(node.wrapped.prev))) :
                                                                                        (attail(node.wrapped) ? nothing : (node.wrapped, HullNode{T}(node.wrapped.next)))
Base.IteratorSize(::HullNodeIterator) = Base.SizeUnknown()


struct PointNodeIterator{T}
    start::PointNode{T}
    rev::Bool
end
"""
    PointNodeIterator(start [, rev])

Returns an iterator over the nodes in the linked list representing the points contained by the convex hull,
starting at the specified node `start`.

If `rev` is `true`, the iterator will advance toward the head of the list.
Otherwise, it will advance toward the tail of the list.
"""
function PointNodeIterator(start::S; rev::Bool = false) where {T,S<:PointNode{T}}
    return PointNodeIterator{T}(start, rev)
end

"""
    PointNodeIterator(hull [, rev])

Returns an iterator over the nodes in the linked list representing the points contained by the convex hull.

If `rev` is `true`, the iterator will start at the tail of the list and advance toward the head.
Otherwise, it will start at the head of the list and advance toward the tail.
"""
PointNodeIterator(h::AbstractConvexHull{T}; rev::Bool = false) where T = PointNodeIterator(rev ? PointNode{T}(h.hull.partner.tail.prev) : PointNode{T}(h.hull.partner.head.next); rev = rev)

Base.iterate(iter::PointNodeIterator) = iterate(iter, iter.start)
Base.iterate(iter::PointNodeIterator{T}, node::S) where {T, S<:PointNode{T}} = iter.rev ? (athead(node.wrapped) ? nothing : (node.wrapped, PointNode{T}(node.wrapped.prev))) :
                                                                                          (attail(node.wrapped) ? nothing : (node.wrapped, PointNode{T}(node.wrapped.next)))
Base.IteratorSize(::PointNodeIterator) = Base.SizeUnknown()


struct BracketedIteratorState{T}
    node::PointNode{T}
    finished::Bool
end

struct BracketedPointNodeIterator{T}
    start::PointNode{T}
    hullstart::HullNode{T}
    hullend::HullNode{T}
    rev::Bool
end
function BracketedPointNodeIterator(start::PointNode{T}, hullstart::HullNode{T}, hullend::HullNode{T}; rev::Bool = false) where T
    hullstart.wrapped.list !== hullend.wrapped.list && throw(ArgumentError("The hull nodes must belong to the same convex hull"))
    start.wrapped.list.partner !== hullstart.wrapped.list && throw(ArgumentError("The starting point node must belong to the same convex hull as the hull points"))
    return BracketedPointNodeIterator{T}(start, hullstart, hullend, rev)
end

Base.iterate(iter::BracketedPointNodeIterator{T}) where T = iterate(iter, BracketedIteratorState{T}(iter.start, iter.start.wrapped.list.len == 0))
function Base.iterate(iter::BracketedPointNodeIterator{T}, state::BracketedIteratorState{T}) where T
    state.finished && return nothing
    if iter.rev
        prevnode = athead(state.node.wrapped.prev) ? tail(state.node.wrapped.list) : state.node.wrapped.prev
        return (state.node.wrapped, BracketedIteratorState(PointNode{T}(prevnode), state.node.wrapped.partner === iter.hullstart.wrapped))
    else
        nextnode = attail(state.node.wrapped.next) ? head(state.node.wrapped.list) : state.node.wrapped.next
        return (state.node.wrapped, BracketedIteratorState(PointNode{T}(nextnode), state.node.wrapped.partner === iter.hullend.wrapped))
    end
end
Base.IteratorSize(::BracketedPointNodeIterator) = Base.SizeUnknown()


function addpoint!(h::AbstractConvexHull{T}, point::T) where T
    # handle the case when the hull is initially empty
    if length(h) == 0
        push!(h.hull, point)
        push!(h.hull.partner, point)
        addpartner!(tail(h.hull), tail(h.hull.partner))
        return h
    end
    if !h.issorted              # if the stored points are unsorted, push the new point to the end
        push!(h.hull.partner, point)
    else                        # otherwise, add as appropriate to maintain sorting
        newpointnode = newnode(h.hull.partner, point)
        f = x -> h.sortedby(x.data) > h.sortedby(point)
        insertbefore = getfirst(f, ListNodeIterator(h.hull.partner))
        isnothing(insertbefore) ? insertnode!(newpointnode, h.hull.partner.tail.prev) : insertnode!(newpointnode, insertbefore.prev)
    end
    if !insidehull(point, h)    # if the new point is outside the hull, update the convex hull
        @show point
        if !h.issorted
            jarvismarch!(h)
        else
            monotonechain!(h)
        end
    end
    return h
end

function addpoints!(h::MutableConvexHull{T}, points::T...; presorted::Bool=false) where T
    h2 = monotonechain(points; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby, presorted=presorted)
    mergehulls(h,h2)
    return h
end

function removepoint!(h::AbstractConvexHull{T}, node::PairedListNode{T}) where T
    (node.list !== h.hull && node.list !== h.hull.partner) && throw(ArgumentError("The specified node must belong to the provided convex hull"))
    partner = node.partner
    deletenode!(node)
    deletenode!(partner)
    length(h) <= 1 && return h
    !h.issorted ? jarvismarch!(h) : monotonechain!(h)
    return h
end



function Base.show(io::IO, h::AbstractConvexHull)
    print(io, typeof(h), '(')
    join(io, h, ", ")
    print(io, ')')
end
