abstract type AbstractChanConvexHull{T} <: AbstractConvexHull{T} end

mutable struct ChanConvexHull{T, F<:Function} <: AbstractChanConvexHull{T}
    hull::TargetedLinkedList{T,PointList{T,HullList{T,F},HullNode{T,HullList{T,F},F},F},PointNode{T,PointList{T,HullList{T,F},HullNode{T,HullList{T,F},F},F},HullNode{T,HullList{T,F},F}}}
    subhulls::Vector{MutableConvexHull{T,F}}
    orientation::HullOrientation
    collinear::Bool
    sortedby::F
end
function ChanConvexHull{T,F}(orientation::HullOrientation=CCW, collinear::Bool=false, sortedby::F=identity) where {T,F}
    subhulls = [MutableConvexHull{T,F}(orientation, collinear, sortedby)]
    hull = TargetedLinkedList{T,PointList{T,HullList{T,F},HullNode{T,HullList{T,F},F},F},PointNode{T,PointList{T,HullList{T,F},HullNode{T,HullList{T,F},F},F},HullNode{T,HullList{T,F},F}}}()
    return ChanConvexHull{T,F}(hull, subhulls, orientation, collinear, sortedby)
end
ChanConvexHull{T}(orientation::HullOrientation=CCW,collinear::Bool=false,sortedby::F=identity) where {T,F} = ChanConvexHull{T,F}(orientation,collinear,sortedby)

function Base.empty!(h::AbstractChanConvexHull) 
    empty!(h.hull)
    subhulls = [subhulls[1]]
    empty!(subhulls[1])
    return h
end

# Iterating a convex hull returns the data contained in the nodes of its hull list
Base.iterate(h::AbstractChanConvexHull) = iterate(h, h.hull.head.next)
Base.iterate(h::AbstractChanConvexHull, node::TargetedListNode) = iterate(h.hull, node)

function addpoint!(h::AbstractChanConvexHull{T}, point::T) where T
    npoints = sum(length, h.subhulls)
    while npoints > 3 && length(h.subhulls)^2 < npoints
        push!(h.subhulls, eltype(h.subhulls)())
    end
    smallhull = argmin(x->length(x),h.subhulls)
    addpoint!(smallhull, point)
    merge_hull_lists!(h)
    return h
end

function mergepoints!(h::AbstractChanConvexHull{T}, points::AbstractVector{T}) where T
    npoints = sum(length, h.subhulls)
    while npoints > 3 && length(h.subhulls)^2 < npoints
        push!(h.subhulls, eltype(h.subhulls)())
    end
    smallhull = argmin(x->length(x),h.subhulls)
    mergepoints!(smallhull, points)
    merge_hull_lists!(h)
    return h
end

function removepoint!(h::AbstractConvexHull{T}, node::TargetedListNode{T}) where T
    node.list !== h.hull && throw(ArgumentError("The specified node must belong to the provided convex hull"))
    shull = getfirst(x -> x === node.target.list, [sh.points for sh in h.subhulls])
    removepoint!(shull, node.target)
    deletenode!(node)
    merge_hull_lists!(h)
    return h
end