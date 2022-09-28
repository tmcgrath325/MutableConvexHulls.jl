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

mutable struct ChanLowerConvexHull{T, F<:Function} <: AbstractChanConvexHull{T}
    hull::TargetedLinkedList{T,PointList{T,HullList{T,F},HullNode{T,HullList{T,F},F},F},PointNode{T,PointList{T,HullList{T,F},HullNode{T,HullList{T,F},F},F},HullNode{T,HullList{T,F},F}}}
    subhulls::Vector{MutableLowerConvexHull{T,F}}
    orientation::HullOrientation
    collinear::Bool
    sortedby::F
end
function ChanLowerConvexHull{T,F}(orientation::HullOrientation=CCW, collinear::Bool=false, sortedby::F=identity) where {T,F}
    subhulls = [MutableLowerConvexHull{T,F}(orientation, collinear, sortedby)]
    hull = TargetedLinkedList{T,PointList{T,HullList{T,F},HullNode{T,HullList{T,F},F},F},PointNode{T,PointList{T,HullList{T,F},HullNode{T,HullList{T,F},F},F},HullNode{T,HullList{T,F},F}}}()
    return ChanLowerConvexHull{T,F}(hull, subhulls, orientation, collinear, sortedby)
end
ChanLowerConvexHull{T}(orientation::HullOrientation=CCW,collinear::Bool=false,sortedby::F=identity) where {T,F} = ChanLowerConvexHull{T,F}(orientation,collinear,sortedby)

mutable struct ChanUpperConvexHull{T, F<:Function} <: AbstractChanConvexHull{T}
    hull::TargetedLinkedList{T,PointList{T,HullList{T,F},HullNode{T,HullList{T,F},F},F},PointNode{T,PointList{T,HullList{T,F},HullNode{T,HullList{T,F},F},F},HullNode{T,HullList{T,F},F}}}
    subhulls::Vector{MutableUpperConvexHull{T,F}}
    orientation::HullOrientation
    collinear::Bool
    sortedby::F
end
function ChanUpperConvexHull{T,F}(orientation::HullOrientation=CCW, collinear::Bool=false, sortedby::F=identity) where {T,F}
    subhulls = [MutableUpperConvexHull{T,F}(orientation, collinear, sortedby)]
    hull = TargetedLinkedList{T,PointList{T,HullList{T,F},HullNode{T,HullList{T,F},F},F},PointNode{T,PointList{T,HullList{T,F},HullNode{T,HullList{T,F},F},F},HullNode{T,HullList{T,F},F}}}()
    return ChanUpperConvexHull{T,F}(hull, subhulls, orientation, collinear, sortedby)
end
ChanUpperConvexHull{T}(orientation::HullOrientation=CCW,collinear::Bool=false,sortedby::F=identity) where {T,F} = ChanUpperConvexHull{T,F}(orientation,collinear,sortedby)

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
        push!(h.subhulls, eltype(h.subhulls)(h.orientation, h.collinear, h.sortedby))
    end
    smallhull = argmin(x->length(x),h.subhulls)
    addpoint!(smallhull, point)
    merge_hull_lists!(h)
    return h
end

function mergepoints!(h::AbstractChanConvexHull{T}, points::AbstractVector{T}) where T
    npoints = sum(length, h.subhulls)
    while npoints > 3 && length(h.subhulls)^2 < npoints
        push!(h.subhulls, eltype(h.subhulls)(h.orientation, h.collinear, h.sortedby))
    end
    smallhull = argmin(x->length(x),h.subhulls)
    mergepoints!(smallhull, points)
    merge_hull_lists!(h)
    return h
end

function removepoint!(h::AbstractChanConvexHull{T}, node::TargetedListNode{T}) where T
    node.list !== h.hull && throw(ArgumentError("The specified node must belong to the provided convex hull"))
    shull = getfirst(x -> x.points === node.target.list, h.subhulls)
    removepoint!(shull, node.target)
    deletenode!(node)
    merge_hull_lists!(h)
    return h
end

function removepoint!(h::AbstractChanConvexHull{T}, node::HullNode{T}) where T
    shull = getfirst(x -> x.hull === node.list, h.subhulls)
    removepoint!(shull, node.target)
    merge_hull_lists!(h)
    return h
end

function removepoint!(h::AbstractChanConvexHull{T}, node::PointNode{T}) where T
    shull = getfirst(x -> x.points === node.list, h.subhulls)
    removepoint!(shull, node.target)
    merge_hull_lists!(h)
    return h
end