abstract type AbstractChanConvexHull{T} <: AbstractConvexHull{T} end

# for debugging
struct ChanHullCache{T}
    data::Vector{T}
    subhulls::Vector{Int}
    removed::Vector{Bool}
end

ChanHullCache{T}() where T = ChanHullCache{T}(T[],Int[],Bool[])

pushcache!(h::H, cache::ChanHullCache{T}, subhull, data::T, removed) where {T, H<:AbstractChanConvexHull{T}} = begin
    shullidx = findfirst(sh -> sh === subhull, h.subhulls)
    push!(cache.data, data)
    push!(cache.subhulls, shullidx)
    push!(cache.removed, removed)
end

pushcache!(h::H, cache::ChanHullCache{T}, subhull, data::AbstractVector{T}, removed) where {T, H<:AbstractChanConvexHull{T}} = begin
    shullidx = findfirst(sh -> sh === subhull, h.subhulls)
    len = length(data)
    append!(cache.data, data)
    append!(cache.subhulls, fill(shullidx, len))
    append!(cache.removed, fill(removed, len))
end

pushcache!(::AbstractChanConvexHull, ::Nothing, subhull, data, removed) = nothing

removecache!(h::AbstractChanConvexHull, cache::ChanHullCache, data, subhull) = begin
    pushcache!(h, cache, subhull, data, true)
end

removecache!(::AbstractChanConvexHull, cache::Nothing, data, subhull) = nothing

emptycache!(cache::ChanHullCache) = begin
    empty!(cache.data)
    empty!(cache.subhulls)
    empty!(cache.removed)
end

emptycache!(::Nothing) = nothing

mutable struct ChanConvexHull{T, F<:Function} <: AbstractChanConvexHull{T}
    hull::TargetedLinkedList{T,PointList{T,HullList{T,F},HullNode{T,HullList{T,F},F},F},PointNode{T,PointList{T,HullList{T,F},HullNode{T,HullList{T,F},F},F},HullNode{T,HullList{T,F},F}}}
    subhulls::Vector{MutableConvexHull{T,F}}
    orientation::HullOrientation
    collinear::Bool
    sortedby::F
    cache::Union{Nothing, ChanHullCache{T}}
    function ChanConvexHull{T,F}(hull, subhulls, orientation, collinear, sortedby, cache=nothing) where {T,F}
        return new(hull, subhulls, orientation, collinear, sortedby, cache)
    end
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
    cache::Union{Nothing, ChanHullCache{T}}
    function ChanLowerConvexHull{T,F}(hull, subhulls, orientation, collinear, sortedby, cache=nothing) where {T,F}
        return new(hull, subhulls, orientation, collinear, sortedby, cache)
    end
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
    cache::Union{Nothing, ChanHullCache{T}}
    function ChanUpperConvexHull{T,F}(hull, subhulls, orientation, collinear, sortedby, cache=nothing) where {T,F}
        new(hull, subhulls, orientation, collinear, sortedby, cache)
    end
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
    emptycache!(h.cache)
    return h
end

# Iterating a convex hull returns the data contained in the nodes of its hull list
Base.iterate(h::AbstractChanConvexHull) = iterate(h, h.hull.head.next)
Base.iterate(h::AbstractChanConvexHull, node::TargetedListNode) = iterate(h.hull, node)

function addpoint!(h::AbstractChanConvexHull{T}, point::T) where T
    npoints = sum(length, h.subhulls)
    while npoints > 3 && length(h.subhulls)^2 < npoints
        push!(h.subhulls, eltype(h.subhulls)(h.orientation, h.collinear, h.sortedby))
        if (h.subhulls[1].points.cache !== nothing)
            h.subhulls[end].points.cache = PairedLinkedLists.SkipListCache{T}()
        end
    end
    smallhull = argmin(x->length(x.points),h.subhulls)
    pushcache!(h, h.cache, smallhull, point, false)
    updatedhull = addpoint!(smallhull, point)[2]
    updatedhull && merge_hull_lists!(h)
    return h, updatedhull
end

function mergepoints!(h::AbstractChanConvexHull{T}, points::AbstractVector{T}) where T
    npoints = sum(length, h.subhulls)
    while npoints > 3 && length(h.subhulls)^2 < npoints
        push!(h.subhulls, eltype(h.subhulls)(h.orientation, h.collinear, h.sortedby))
        if (h.subhulls[1].points.cache !== nothing)
            h.subhulls[end].points.cache = PairedLinkedLists.SkipListCache{T}()
        end
    end
    smallhull = argmin(x->length(x.points),h.subhulls)
    pushcache!(h, h.cache, smallhull, points, false)
    mergepoints!(smallhull, points)
    merge_hull_lists!(h)
    return h
end

function removepoint!(h::AbstractChanConvexHull{T}, node::TargetedListNode{T}) where T
    updatedhull = removepoint!(h, node.target)[2]
    return h, updatedhull
end

function removepoint!(h::AbstractChanConvexHull{T}, node::HullNode{T}) where T
    updatedhull = removepoint!(h, node.target)[2]
    return h, updatedhull
end

function removepoint!(h::AbstractChanConvexHull{T}, node::PointNode{T}) where T
    shull = getfirst(sh -> sh.points === node.list, h.subhulls)
    shull === nothing && throw(ArgumentError("The specified node must belong to the provided convex hull"))
    updatedhull = removepoint!(shull, node.target)[2]
    updatedhull && merge_hull_lists!(h)
    removecache!(h, h.cache, node.data, shull)
    return h, updatedhull
end


function copyfromcache(h::H) where {T,H<:AbstractChanConvexHull{T}}
    @assert !isnothing(h.cache)
    hcopy = H(h.orientation, h.collinear, h.sortedby)
    hcopy.cache = ChanHullCache{T}()
    hcopy.subhulls[1].points.cache = typeof(h.subhulls[1].points.cache)()
    for i=2:length(h.subhulls)
        push!(hcopy.subhulls, eltype(h.subhulls)(h.orientation, h.collinear, h.sortedby))
        hcopy.subhulls[end].points.cache = typeof(h.subhulls[i].points.cache)()
    end
    pointstoadd = T[]
    levelstoadd = Int[]
    shullcounters = fill(0, length(h.subhulls))
    currentshullidx = first(h.cache.subhulls)
    for (i,(point, shullidx, removed)) in enumerate(zip(h.cache.data, h.cache.subhulls, h.cache.removed))
        if !isempty(pointstoadd) && (shullidx != currentshullidx || removed || i == length(h.cache.data))
            if isempty(levelstoadd)
                mergepoints!(hcopy.subhulls[currentshullidx], pointstoadd)
            else
                pointstoadd = sort(pointstoadd; by = h.sortedby)
                @assert pointstoadd == h.subhulls[currentshullidx].points.cache.data[(length(hcopy.subhulls[currentshullidx].points.cache.data)+1):shullcounters[currentshullidx]]
                @assert levelstoadd == h.subhulls[currentshullidx].points.cache.levels[(length(hcopy.subhulls[currentshullidx].points.cache.levels)+1):shullcounters[currentshullidx]]
                for (point, level) in zip(pointstoadd, levelstoadd)
                    PairedLinkedLists.pushskip!(hcopy.subhulls[currentshullidx].points, point, level)
                    monotonechain!(hcopy.subhulls[currentshullidx])
                end
            end
            pushcache!(hcopy, hcopy.cache, hcopy.subhulls[currentshullidx], pointstoadd, false)
            merge_hull_lists!(hcopy)
            empty!(pointstoadd)
            empty!(levelstoadd)
        end
        currentshullidx = shullidx
        if removed
            node = getfirst(x -> x.data == point, PointNodeIterator(hcopy.subhulls[shullidx].points.head.next))
            removepoint!(hcopy, node)
            shullcounters[shullidx] += 1
        else
            push!(pointstoadd, point)
            if (h.subhulls[shullidx].points.cache !== nothing)
                push!(levelstoadd, h.subhulls[shullidx].points.cache.levels[shullcounters[shullidx]+1])
            end
            shullcounters[shullidx] += 1
        end
    end
    return hcopy
end


function chanhullsidentical(h1::AbstractChanConvexHull, h2::AbstractChanConvexHull)
    h1.hull == h2. hull || return false
    for (sh1, sh2) in zip(h1.subhulls, h2.subhulls)
        sh1 == sh2 || return false
        sh1.points == sh2.points || return false
        PairedLinkedLists.skiplistsidentical(sh1.points, sh2.points) || return false
    end
    return true
end