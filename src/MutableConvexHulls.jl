"""
MutableConvexHulls.jl provides ways to calculate and update vertex representations of planar [convex polytopes](https://en.wikipedia.org/wiki/Convex_polytope) 
(i.e. convex hulls). It is intended for use in situations when a convex hull must be updated iteratively with addition or removal of points.

See also [MutableConvexHull](@ref), [MutableLowerConvexHull](@ref), [MutableUpperConvexHull](@ref), [monotonechain](@ref), [addpoint!](@ref), [mergepoints!](@ref), [removepoint!](@ref)
"""
module MutableConvexHulls

using DoubleFloats: DoubleFloats, DoubleFloat
using PairedLinkedLists: PairedLinkedLists,
    AbstractNode, AbstractList,
    AbstractPairedListNode, AbstractPairedSkipNode,
    AbstractPairedLinkedList, AbstractPairedSkipList,
    SkipListCache,
    ListNodeIterator,
    TargetedLinkedList, TargetedListNode,
    addtarget!, athead, attail, deletenode!, hastarget,
    head, insertafter!, newnode, nodetype, removetarget!, search, tail

include("utils.jl")
include("orientation.jl")
include("pointlist.jl")
include("convexhull.jl")
include("chanhull.jl")
include("monotonechain.jl")
include("jarvismarch.jl")
include("inside.jl")
include("merge.jl")

export AbstractConvexHull, MutableConvexHull, MutableLowerConvexHull, MutableUpperConvexHull
export AbstractChanConvexHull, ChanConvexHull, ChanLowerConvexHull, ChanUpperConvexHull
export HullList, PointList, HullNode, PointNode
export addpoint!, mergepoints!, removepoint!
export monotonechain, lower_monotonechain, upper_monotonechain
export jarvismarch, lower_jarvismarch, upper_jarvismarch
export CCW, CW
export insidehull
export mergehulls, mergehulls!

# `HullNodeIterator` and `PointNodeIterator` expose `HullNode`/`PointNode` field
# layouts, so they are public rather than exported. The `public` keyword exists
# only on Julia ≥ 1.11; on the LTS these remain reachable as qualified names.
@static if VERSION >= v"1.11"
    eval(Expr(:public, :HullNodeIterator, :PointNodeIterator))
end

end
