# cross product of 2D vectors OA and OB
cross2d(o, a, b) = cross2d(triplet_to_vectors(o, a, b)...)
cross2d(vec1, vec2) = (vec1[1] * vec2[2] - vec1[2] * vec2[1])

dot2d(vec1, vec2) = (vec1[1] * vec2[1] + vec1[2] * vec2[2])