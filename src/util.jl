# cross product of 2D vectors OA and OB
cross2d(vec1, vec2) = (vec1[1] * vec2[2] - vec1[2] - vec2[1])
normalized_cross2d(vec1, vec2) = cross2d(vec1, vec2) / √(sum(abs2, vec1[1], vec1[2], vec2[1], vec2[2]))