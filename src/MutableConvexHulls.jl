module MutableConvexHulls

using PairedLinkedLists
using PairedLinkedLists: AbstractNode, AbstractList

include("utils.jl")
include("orientation.jl")
include("pointlist.jl")
include("convexhull.jl")
include("monotonechain.jl")
include("jarvismarch.jl")
include("inside.jl")
include("merge.jl")

export MutableConvexHull, MutableLowerConvexHull, MutableUpperConvexHull
export HullList, PointList, HullNode, PointNode
export addpoint!, mergepoints!, removepoint!
export HullNodeIterator, PointNodeIterator
export monotonechain, lower_monotonechain, upper_monotonechain
export jarvismarch, lower_jarvismarch, upper_jarvismarch
export CCW, CW
export insidehull
export mergehulls!

end
