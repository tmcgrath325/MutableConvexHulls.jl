"""
    AbstractConvexHull{T}

Abstract supertype for all convex hull types over 2-D points of type `T`.

Concrete subtypes: [`MutableConvexHull`](@ref), [`MutableLowerConvexHull`](@ref),
[`MutableUpperConvexHull`](@ref), [`ChanConvexHull`](@ref), [`ChanLowerConvexHull`](@ref),
[`ChanUpperConvexHull`](@ref).
"""
abstract type AbstractConvexHull{T} end

mutable struct MutableConvexHull{T, F <: Function} <: AbstractConvexHull{T}
    const hull::HullList{T, F}
    const points::PointList{T, HullList{T, F}, HullNode{T, HullList{T, F}, F}, F}
    const orientation::HullOrientation
    const collinear::Bool
    const sortedby::F
end
function MutableConvexHull{T, F}(; orientation::HullOrientation = CCW, collinear::Bool = false, sortedby::F = identity) where {T, F}
    points = PointList{T}(; sortedby = sortedby)
    hull = HullList{T, F}()
    addtarget!(hull, points)
    return MutableConvexHull{T, F}(hull, points, orientation, collinear, sortedby)
end

"""
    h = MutableConvexHull{T}(; orientation=CCW, collinear=false, sortedby=identity)

A mutable convex hull over 2-D points of type `T` (e.g., `Tuple{Float64,Float64}`),
supporting incremental addition and removal of points. Iterating `h` yields the
hull-vertex sequence in the specified orientation.

Initialize an empty hull with the provided attributes.

`orientation` specifies whether the points along the convex hull are ordered clockwise `CW`, or counterclockwise `CCW`, and defaults to `CCW`.

`collinear` specifies whether collinear points are allowed along the surface of the convex hull, and defaults to `false`.

`sortedby` specifies a function to apply to points prior to sorting, and defaults to `identity` (resulting in default sorting behavior).

See also: [monotonechain](@ref), [jarvismarch](@ref), [addpoint!](@ref), [mergepoints!](@ref), [removepoint!](@ref)
"""
MutableConvexHull{T}(; orientation = CCW, collinear = false, sortedby::F = identity) where {T, F} = MutableConvexHull{T, F}(; orientation, collinear, sortedby)

mutable struct MutableLowerConvexHull{T, F <: Function} <: AbstractConvexHull{T}
    const hull::HullList{T, F}
    const points::PointList{T, HullList{T, F}, HullNode{T, HullList{T, F}, F}, F}
    const orientation::HullOrientation
    const collinear::Bool
    const sortedby::F
end
function MutableLowerConvexHull{T, F}(; orientation::HullOrientation = CCW, collinear::Bool = false, sortedby::F = identity) where {T, F}
    points = PointList{T}(; sortedby = sortedby)
    hull = HullList{T, F}()
    addtarget!(hull, points)
    return MutableLowerConvexHull{T, F}(hull, points, orientation, collinear, sortedby)
end

"""
    h = MutableLowerConvexHull{T}(; orientation=CCW, collinear=false, sortedby=identity)

A mutable lower convex hull over 2-D points of type `T` (e.g., `Tuple{Float64,Float64}`):
the chain of hull vertices from the leftmost to the rightmost point along the bottom boundary.
Supports the same incremental operations as [`MutableConvexHull`](@ref).

Initialize an empty hull with the provided attributes.

`orientation` specifies whether the points along the convex hull are ordered clockwise `CW`, or counterclockwise `CCW`, and defaults to `CCW`.

`collinear` specifies whether collinear points are allowed along the surface of the convex hull, and defaults to `false`.

`sortedby` specifies a function to apply to points prior to sorting, and defaults to `identity` (resulting in default sorting behavior).

See also: [lower_monotonechain](@ref), [lower_jarvismarch](@ref), [addpoint!](@ref), [mergepoints!](@ref), [removepoint!](@ref)
"""
MutableLowerConvexHull{T}(; orientation = CCW, collinear = false, sortedby::F = identity) where {T, F} = MutableLowerConvexHull{T, F}(; orientation, collinear, sortedby)

mutable struct MutableUpperConvexHull{T, F <: Function} <: AbstractConvexHull{T}
    const hull::HullList{T, F}
    const points::PointList{T, HullList{T, F}, HullNode{T, HullList{T, F}, F}, F}
    const orientation::HullOrientation
    const collinear::Bool
    const sortedby::F
end
function MutableUpperConvexHull{T, F}(; orientation::HullOrientation = CCW, collinear::Bool = false, sortedby::F = identity) where {T, F}
    points = PointList{T}(; sortedby = sortedby)
    hull = HullList{T, F}()
    addtarget!(hull, points)
    return MutableUpperConvexHull{T, F}(hull, points, orientation, collinear, sortedby)
end

"""
    h = MutableUpperConvexHull{T}(; orientation=CCW, collinear=false, sortedby=identity)

A mutable upper convex hull over 2-D points of type `T` (e.g., `Tuple{Float64,Float64}`):
the chain of hull vertices from the leftmost to the rightmost point along the top boundary.
Supports the same incremental operations as [`MutableConvexHull`](@ref).

Initialize an empty hull with the provided attributes.

`orientation` specifies whether the points along the convex hull are ordered clockwise `CW`, or counterclockwise `CCW`, and defaults to `CCW`.

`collinear` specifies whether collinear points are allowed along the surface of the convex hull, and defaults to `false`.

`sortedby` specifies a function to apply to points prior to sorting, and defaults to `identity` (resulting in default sorting behavior).

See also: [upper_monotonechain](@ref), [upper_jarvismarch](@ref), [addpoint!](@ref), [mergepoints!](@ref), [removepoint!](@ref)
"""
MutableUpperConvexHull{T}(; orientation = CCW, collinear = false, sortedby::F = identity) where {T, F} = MutableUpperConvexHull{T, F}(; orientation, collinear, sortedby)

Base.isempty(h::AbstractConvexHull) = isempty(h.hull)
Base.length(h::AbstractConvexHull) = length(h.hull)
Base.eltype(h::AbstractConvexHull) = eltype(h.hull)

Base.:(==)(h1::AbstractConvexHull, h2::AbstractConvexHull) = h1.hull == h2.hull

# Hash the hull-vertex sequence that `==` compares, so equal hulls hash equally
# and behave correctly as `Dict` keys or `Set` elements.
function Base.hash(h::AbstractConvexHull, x::UInt)
    x = hash(:AbstractConvexHull, x)
    for data in h
        x = hash(data, x)
    end
    return x
end

Base.empty!(h::AbstractConvexHull) = empty!(h.hull)
Base.empty(h::H) where {H <: AbstractConvexHull} = H(; orientation = h.orientation, collinear = h.collinear, sortedby = h.sortedby)

# Deep copy by replaying the contained points: the result shares no linked-list
# nodes with `h`, so mutating one hull never affects the other.
function Base.copy(h::H) where {H <: AbstractConvexHull}
    hcopy = H(; orientation = h.orientation, collinear = h.collinear, sortedby = h.sortedby)
    for node in PointNodeIterator(h)
        addpoint!(hcopy, node.data)
    end
    return hcopy
end

# Iterating a convex hull returns the data contained in the nodes of its hull list
Base.iterate(h::AbstractConvexHull) = iterate(h, h.hull.head.next)
Base.iterate(h::AbstractConvexHull, node::HullNode) = iterate(h.hull, node)

struct HullNodeIterator{S <: HullNode}
    start::S
    rev::Bool
end
"""
    HullNodeIterator(start[; rev=false])

Return an iterator over the [`HullNode`](@ref) elements in the convex hull's linked
list, starting at `start`.

If `rev` is `true`, the iterator advances toward the head of the list.
Otherwise, it advances toward the tail.
"""
function HullNodeIterator(start::S; rev::Bool = false) where {S <: HullNode}
    return HullNodeIterator{S}(start, rev)
end

"""
    HullNodeIterator(hull[; rev=false])

Return an iterator over the [`HullNode`](@ref) elements in the convex hull's linked list.

If `rev` is `true`, iteration starts at the tail and advances toward the head.
Otherwise, it starts at the head and advances toward the tail.
"""
HullNodeIterator(h::AbstractConvexHull; rev::Bool = false) = HullNodeIterator{nodetype(h.hull)}(rev ? h.hull.tail.prev : h.hull.head.next, rev)

Base.iterate(iter::HullNodeIterator) = iterate(iter, iter.start)
Base.iterate(iter::HullNodeIterator{S}, node::S) where {S <: HullNode} = iter.rev ? (athead(node) ? nothing : (node, node.prev)) :
    (attail(node) ? nothing : (node, node.next))
Base.IteratorSize(::HullNodeIterator) = Base.SizeUnknown()
Base.IteratorEltype(::Type{<:HullNodeIterator}) = Base.HasEltype()
Base.eltype(::Type{<:HullNodeIterator{S}}) where {S <: HullNode} = S


struct PointNodeIterator{S <: PointNode}
    start::S
    rev::Bool
end
"""
    PointNodeIterator(start[; rev=false])

Return an iterator over the [`PointNode`](@ref) elements in the hull's point list,
starting at `start`.

If `rev` is `true`, the iterator advances toward the head of the list.
Otherwise, it advances toward the tail.
"""
function PointNodeIterator(start::S; rev::Bool = false) where {S}
    return PointNodeIterator{S}(start, rev)
end

"""
    PointNodeIterator(hull[; rev=false])

Return an iterator over the [`PointNode`](@ref) elements in the hull's point list.

If `rev` is `true`, iteration starts at the tail and advances toward the head.
Otherwise, it starts at the head and advances toward the tail.
"""
PointNodeIterator(h::AbstractConvexHull; rev::Bool = false) = PointNodeIterator{nodetype(h.points)}(rev ? h.points.tail.prev : h.points.head.next, rev)

Base.iterate(iter::PointNodeIterator) = iterate(iter, iter.start)
Base.iterate(iter::PointNodeIterator{S}, node::S) where {S} = iter.rev ? (athead(node) ? nothing : (node, node.prev)) :
    (attail(node) ? nothing : (node, node.next))
Base.IteratorSize(::PointNodeIterator) = Base.SizeUnknown()
Base.IteratorEltype(::Type{<:PointNodeIterator}) = Base.HasEltype()
Base.eltype(::Type{<:PointNodeIterator{S}}) where {S <: PointNode} = S

"""
    addpoint!(hull, point)

Add `point` to the list of points contained by the provided convex hull `hull`. If `point` lies outside the convex hull,
the list of hull points will be updated accordingly.

Return `(hull, expanded)` where `expanded` is `true` if the added point expands the hull
and `false` if it lies inside or on the boundary.

# Examples
```jldoctest
julia> h = MutableConvexHull{Tuple{Float64,Float64}}();

julia> h, _ = addpoint!(h, (0.0, 0.0)); h, _ = addpoint!(h, (1.0, 0.0));

julia> h, expanded = addpoint!(h, (0.0, 1.0));

julia> expanded
true

julia> h, expanded = addpoint!(h, (0.25, 0.25));

julia> expanded
false
```

See also: [mergepoints!](@ref), [removepoint!](@ref)
"""
function addpoint!(h::AbstractConvexHull{T}, point::T) where {T}
    # handle the case when the hull is initially empty
    if length(h) == 0
        push!(h.hull, point)
        push!(h.points, point)
        addtarget!(tail(h.hull), tail(h.points))
        return h, true
    end
    push!(h.points, point)
    if !insidehull(point, h)    # if the new point is outside the hull, update the convex hull
        monotonechain!(h)
        return h, true
    end
    return h, false
end

"""
    mergepoints!(hull, points)

Add `points` to the list of points contained by the provided convex hull `hull`. If any of the `points` lie outside the convex hull,
the list of hull points will be updated accordingly.

`points` may be a vector of points or an `AbstractMatrix` in which each row is one point.

This function finds the convex hull of the `points` to be added before merging them with the `hull`. See
[Chan's algorithm](https://en.wikipedia.org/wiki/Chan%27s_algorithm) for a similar idea.

Return `hull`.

See also: [addpoint!](@ref), [removepoint!](@ref)
"""
function mergepoints!(h::MutableConvexHull{T}, points::AbstractVector{T}) where {T}
    h2 = monotonechain(points; orientation = h.orientation, collinear = h.collinear, sortedby = h.sortedby)
    mergehulls!(h, h2)
    return h
end
function mergepoints!(h::MutableLowerConvexHull{T}, points::AbstractVector{T}) where {T}
    h2 = lower_monotonechain(points; orientation = h.orientation, collinear = h.collinear, sortedby = h.sortedby)
    mergehulls!(h, h2)
    return h
end
function mergepoints!(h::MutableUpperConvexHull{T}, points::AbstractVector{T}) where {T}
    h2 = upper_monotonechain(points; orientation = h.orientation, collinear = h.collinear, sortedby = h.sortedby)
    mergehulls!(h, h2)
    return h
end
mergepoints!(h::AbstractConvexHull, points::AbstractMatrix) = mergepoints!(h, rowpoints(points))

"""
    removepoint!(hull, node)

Remove `node` from `hull`. `node` may be a [`HullNode`](@ref) (a vertex on the hull
boundary) or a [`PointNode`](@ref) (any tracked point, inside or on the boundary). If
removing `node` alters the hull boundary, it is recomputed.

Return `(hull, removed)` where `removed` is `true` if a hull-boundary vertex was
removed and `false` if an interior or duplicate point was removed.

See also: [addpoint!](@ref), [mergepoints!](@ref)
"""
function removepoint!(h::AbstractConvexHull{T}, node::HullNode{T}) where {T}
    updatedhull = removepoint!(h, node.target)[2]
    return h, updatedhull
end

function removepoint!(h::Union{MutableLowerConvexHull{T}, MutableUpperConvexHull{T}}, node::PointNode{T}) where {T}
    node.list !== h.points && throw(ArgumentError("The specified node must belong to the provided convex hull"))
    if hastarget(node)
        # handle cases with duplicate data
        next = node.next
        prev = node.prev
        if !athead(prev)
            if coordsareequal(prev.data, node.data)
                removepoint!(h, prev)
                return h, false
            end
        end
        if !attail(next)
            if coordsareequal(next.data, node.data)
                removepoint!(h, next)
                return h, false
            end
        end
        # handle general case
        target = node.target
        hullstart = target.prev
        hullstop = target.next
        deletenode!(target)
        deletenode!(node)
        start = athead(hullstart) ? firstpoint(h) : hullstart.target
        stop = attail(hullstop) ? lastpoint(h) : hullstop.target
        coordsareequal(start.data, stop.data) && return (h, true)
        monotonechain!(h, start, stop)
        return h, true
    else
        deletenode!(node)
        return h, false
    end
end

function removepoint!(h::MutableConvexHull{T}, node::PointNode{T}) where {T}
    node.list !== h.points && throw(ArgumentError("The specified node must belong to the provided convex hull"))
    if hastarget(node)
        # handle cases with duplicate data
        next = node.next
        prev = node.prev
        if !athead(prev)
            if coordsareequal(prev.data, node.data)
                removepoint!(h, prev)
                return h, false
            end
        end
        if !attail(next)
            if coordsareequal(next.data, node.data)
                removepoint!(h, next)
                return h, false
            end
        end
        # handle general case
        target = node.target
        deletenode!(target)
        deletenode!(node)
        monotonechain!(h)
        return h, true
    else
        deletenode!(node)
        return h, false
    end
end

# Locate the point node whose data equals `value`, or `nothing` if absent.
# `search` lands on the last node whose sort key is ≤ that of `value`; scanning
# back across the equal-key run and comparing full coordinates handles a
# `sortedby` that maps distinct points to the same key.
function findpointnode(points::PointList{T}, value::T) where {T}
    node = search(points, value)
    svalue = points.sortedby(value)
    while !athead(node) && points.sortedby(node.data) == svalue
        coordsareequal(node.data, value) && return node
        node = node.prev
    end
    return nothing
end

"""
    removepoint!(hull, value)

Locate the point equal to `value` in `hull` and remove it, delegating to
[`removepoint!`](@ref)`(hull, node)`. Throws an `ArgumentError` if no point
equal to `value` is contained in `hull`.

Return `(hull, removed)` where `removed` is `true` if a hull-boundary vertex
was removed and `false` if an interior or duplicate point was removed.

# Examples
```jldoctest
julia> h = monotonechain([(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (0.25, 0.25)]);

julia> h, removed = removepoint!(h, (0.25, 0.25));

julia> removed  # interior point; hull boundary unchanged
false

julia> h, removed = removepoint!(h, (1.0, 0.0));

julia> removed  # hull vertex; boundary recomputed
true

julia> h
MutableConvexHull{Tuple{Float64, Float64}, typeof(identity)}((0.0, 0.0), (0.0, 1.0))
```

See also: [addpoint!](@ref), [mergepoints!](@ref)
"""
function removepoint!(h::AbstractConvexHull{T}, value::T) where {T}
    node = findpointnode(h.points, value)
    node === nothing && throw(ArgumentError("No point equal to $value is contained in the convex hull"))
    return removepoint!(h, node)
end

function Base.show(io::IO, h::AbstractConvexHull)
    print(io, typeof(h), '(')
    join(io, h, ", ")
    return print(io, ')')
end
