abstract type AbstractConvexHull{T} end

struct MutableConvexHull{T} <: AbstractConvexHull{T}
    hull::PairedLinkedList{T}
    orientation::HullOrientation
    collinear::Bool
end
function MutableConvexHull{T}(orientation::HullOrientation = CCW, collinear::Bool = false) where T
    pointslist = PairedLinkedList{T}()
    hull = PairedLinkedList{T}()
    addpartner!(hull,pointslist)
    return MutableConvexHull{T}(hull, orientation, collinear)
end

struct MutableLowerConvexHull{T} <: AbstractConvexHull{T}
    hull::PairedLinkedList{T}
    orientation::HullOrientation
    collinear::Bool
end
function MutableLowerConvexHull{T}(orientation::HullOrientation = CCW, collinear::Bool = false) where T
    pointslist = PairedLinkedList{T}()
    hull = PairedLinkedList{T}()
    addpartner!(hull,pointslist)
    return MutableLowerConvexHull{T}(hull, orientation, collinear)
end

struct MutableUpperConvexHull{T} <: AbstractConvexHull{T}
    hull::PairedLinkedList{T}
    orientation::HullOrientation
    collinear::Bool
end
function MutableUpperConvexHull{T}(orientation::HullOrientation = CCW, collinear::Bool = false) where T
    pointslist = PairedLinkedList{T}()
    hull = PairedLinkedList{T}()
    addpartner!(hull,pointslist)
    return MutableUpperConvexHull{T}(hull, orientation, collinear)
end

Base.length(h::AbstractConvexHull) = length(h.hull)

function Base.copy(h::H) where {T,H<:AbstractConvexHull{T}}
    copiedhull = copy(h.hull)
    return H(copiedhull, h.orientation, h.collinear)
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

