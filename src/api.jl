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