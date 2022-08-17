abstract type AbstractConvexHull{T} end

struct MutableConvexHull{T} <: AbstractConvexHull{T}
    hull::PairedLinkedList{T}
    orientation::HullOrientation
    collinear::Bool
    sortedby::Union{Nothing, Function}
end
function MutableConvexHull{T}(orientation::HullOrientation = CCW, collinear::Bool = false, sortedby::Union{Nothing,Function}=nothing) where T
    pointslist = PairedLinkedList{T}()
    hull = PairedLinkedList{T}()
    addpartner!(hull,pointslist)
    return MutableConvexHull{T}(hull, orientation, collinear, sortedby)
end

struct MutableLowerConvexHull{T} <: AbstractConvexHull{T}
    hull::PairedLinkedList{T}
    orientation::HullOrientation
    collinear::Bool
    sortedby::Union{Nothing, Function}
end
function MutableLowerConvexHull{T}(orientation::HullOrientation = CCW, collinear::Bool = false, sortedby::Union{Nothing,Function}=nothing) where T
    pointslist = PairedLinkedList{T}()
    hull = PairedLinkedList{T}()
    addpartner!(hull,pointslist)
    return MutableLowerConvexHull{T}(hull, orientation, collinear, sortedby)
end

struct MutableUpperConvexHull{T} <: AbstractConvexHull{T}
    hull::PairedLinkedList{T}
    orientation::HullOrientation
    collinear::Bool
    sortedby::Union{Nothing, Function}
end
function MutableUpperConvexHull{T}(orientation::HullOrientation = CCW, collinear::Bool = false) where T
    pointslist = PairedLinkedList{T}()
    hull = PairedLinkedList{T}()
    addpartner!(hull,pointslist)
    return MutableUpperConvexHull{T}(hull, orientation, collinear)
end

Base.length(h::AbstractConvexHull) = length(h.hull)

function Base.copy(h::H) where {T,H<:AbstractConvexHull{T}}
    hcopy = H(h.orientation, h.collinear)
    pointslist = hcopy.hull.partner
    hulltargets = TargetedLinkedList(pointslist)
    for pointsnode in IteratingListNodes(h.hull.partner)
        push!(pointslist, pointsnode.data)
        if haspartner(pointsnode) 
            push!(hulltargets, pointsnode.data)
            addpartner!(hulltargets.tail.prev, pointslist.tail.prev)
        end
    end
    for hullnode in IteratingListNodes(h.hull)
        push!(hcopy.hull, hullnode.data)
        targetingnode = first(Iterators.filter(x->x.data == hullnode.data, IteratingListNodes(hulltargets)))
        addpartner!(hcopy.hull.tail.prev, targetingnode.partner)
    end
    return hcopy
end

struct PointNode{T<:AbstractListNode}
    wrapped::T
end
struct HullNode{T<:AbstractListNode}
    wrapped::T
end

struct IteratingPointsListNodes{S<:AbstractListNode}
    start::S
    rev::Bool
    function IteratingPointsListNodes(hull::AbstractConvexHull; rev::Bool = false)
        return new{PairedLinkedLists.nodetype(hull.hull)}(hull.hull.partner.head.next, rev)
    end
end

struct IteratingHullListNodes{S<:AbstractListNode}
    start::S
    rev::Bool
    function IteratingHullListNodes(hull::AbstractConvexHull; rev::Bool = false)
        return new{PairedLinkedLists.nodetype(hull.hull)}(hull.hull.head.next, rev)
    end
end

Base.iterate(h::AbstractConvexHull) = Base.iterate(h.hull)

