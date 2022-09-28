"""
MutableConvexHulls.jl provides ways to calculate and update vertex representations of planar [convex polytopes](https://en.wikipedia.org/wiki/Convex_polytope) 
(i.e. convex hulls). It is intended for use in situations when a convex hull must be updated iteratively with addition or removal of points.

See also [MutableConvexHull](@ref), [LowerMutableConvexHull](@ref), [UpperMutableConvexHull](@ref), [monotonechain](@ref), [addpoint!](@ref), [mergepoints!](@ref), [removepoint!](@ref), 
"""
module MutableConvexHulls

using PairedLinkedLists
using PairedLinkedLists: AbstractNode, AbstractList

include("utils.jl")
include("orientation.jl")
include("pointlist.jl")
include("convexhull.jl")
include("chanhull.jl")
include("monotonechain.jl")
include("jarvismarch.jl")
include("inside.jl")
include("merge.jl")

export MutableConvexHull, MutableLowerConvexHull, MutableUpperConvexHull
export ChanConvexHull, ChanLowerConvexHull, ChanUpperConvexHull
export HullList, PointList, HullNode, PointNode
export addpoint!, mergepoints!, removepoint!
export HullNodeIterator, PointNodeIterator
export monotonechain, lower_monotonechain, upper_monotonechain
export jarvismarch, lower_jarvismarch, upper_jarvismarch
export CCW, CW
export insidehull
export mergehulls!

end
