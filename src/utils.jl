# cross product of 2D vectors OA and OB
# cross2d(o, a, b) = cross2d(a .- o, b .- o)
cross2d(oa, ob) = (oa[1] * ob[2] - oa[2] * ob[1])

# dot2d(o, a, b) = dot2d(a .- o, b .- o)
dot2d(oa, ob) = (oa[1] * ob[1] + oa[2] * ob[2])