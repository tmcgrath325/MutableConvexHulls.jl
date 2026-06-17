# planar coordinates
coordsareequal(data1, data2) = data1[1] == data2[1] && data1[2] == data2[2]

# Treat each row of a matrix as one point, returning a vector of point tuples.
# Iterates `eachrow` so rows are taken in `axes(points, 1)` order, supporting
# views, adjoints, and matrices with non-1-based axes.
rowpoints(points::AbstractMatrix) = [(row...,) for row in eachrow(points)]

# 2D subtraction of data
sub2d(a, b) = (DoubleFloat(a[1]) - DoubleFloat(b[1]), DoubleFloat(a[2]) - DoubleFloat(b[2]))

# cross and dot products of 2D vectors OA and OB
cross2d(oa, ob) = (oa[1] * ob[2] - oa[2] * ob[1])
dot2d(oa, ob) = (oa[1] * ob[1] + oa[2] * ob[2])

# linear interpolation
function linterp(x, start, finish)
    xratio = (x - start[1]) / (finish[1] - start[1])
    y = start[2] + xratio * (finish[2] - start[2])
    return y
end

# returns the first element of the iterator that satisfies the predicate
function getfirst(p, itr)
    for el in itr
        p(el) && return el
    end
    return nothing
end