using MutableConvexHulls
using Random

Random.seed!(1234)

orientation = CW
collinear = true
sortedby = x -> (x[1], -x[2])
coords = [(1, 5), (4, 2)]

h = MutableConvexHull{eltype(coords)}(orientation, collinear, sortedby)
for (i,coord) in enumerate(coords[1:end-1])
    addpoint!(h, coord)
    @show coord
    @assert h == jarvismarch(coords[1:i]; orientation=orientation, collinear=collinear, sortedby=sortedby)
end

@show jarvismarch(coords; orientation=orientation, collinear=collinear, sortedby=sortedby)
addpoint!(h, coords[end])