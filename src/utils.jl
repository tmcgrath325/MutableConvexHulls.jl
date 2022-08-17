# cross product of 2D vectors OA and OB
# cross2d(o, a, b) = cross2d(a .- o, b .- o)
cross2d(oa, ob) = (oa[1] * ob[2] - oa[2] * ob[1])

# dot2d(o, a, b) = dot2d(a .- o, b .- o)
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