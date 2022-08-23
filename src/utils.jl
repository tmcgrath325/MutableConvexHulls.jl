# check planar coordinates for equality
coordsareequal(data1, data2) = data1[1] == data2[1] && data1[2] == data2[2]

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