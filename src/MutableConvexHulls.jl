module MutableConvexHulls

using PairedLinkedLists

include("utils.jl")
include("orientation.jl")
include("api.jl")
include("monotonechain.jl")
include("jarvismarch.jl")
include("inside.jl")
include("chan.jl")

export MutableConvexHull, MutableLowerConvexHull, MutableUpperConvexHull
export addpoint!, removepoint!
export HullNode, PointNode, HullNodeIterator, PointNodeIterator, BracketedPointNodeIterator
export monotonechain, lower_monotonechain, upper_monotonechain
export jarvismarch, lower_jarvismarch, upper_jarvismarch
export CCW, CW
export insidehull
export mergehulls!, mergehulls

end
