abstract type AbstractConvexHull{T} end

struct MutableConvexHull{T} <: AbstractConvexHull{T}
    hull::PairedLinkedList{T}
    orientation::HullOrientation
    colinear::Bool
end

struct MutableLowerConvexHull{T} <: AbstractConvexHull{T}
    hull::PairedLinkedList{T}
    orientation::HullOrientation
    colinear::Bool
end

struct MutableUpperConvexHull{T} <: AbstractConvexHull{T}
    hull::PairedLinkedList{T}
    orientation::HullOrientation
    colinear::Bool
end

Base.length(h::AbstractConvexHull) = length(h.hull)

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
