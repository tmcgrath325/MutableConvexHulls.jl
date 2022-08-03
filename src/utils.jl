# clockwise/counterclockwise selection
abstract type HullDirection end
struct CW <: HullDirection end
struct CCW <: HullDirection end

# from three points O, A, and B, return vectors OA and OB
triplet_to_vectors(o, a, b) = (a.-o, b.-o)

# cross product of 2D vectors OA and OB
cross2d(o, a, b) = cross2d(triplet_to_vectors(o, a, b)...)
cross2d(vec1, vec2) = (vec1[1] * vec2[2] - vec1[2] * vec2[1])

dot2d(vec1, vec2) = (vec1[1] * vec2[1] + vec1[2] * vec2[2])

# determine if the change from one vector to another represents a left (ccw) or right (cw) turn. 
# For the methods below, colinear vectors will always yield `true` in both cases.
left_turn(o, a, b) = left_turn(triplet_to_vectors(o, a, b)...)
left_turn(vec1, vec2) = cross2d(vec1, vec2) >= 0
right_turn(o, a, b) = right_turn(triplet_to_vectors(o, a, b)...)
right_turn(vec1, vec2) = cross2d(vec1, vec2) <= 0 

# For the methods below, colinear vectors will only yield `true` if they are aligned (dot product > 0)
aligned_left_turn(o, a, b) = aligned_left_turn(triplet_to_vectors(o, a, b)...)
function aligned_left_turn(vec1, vec2)
    cp = cross2d(vec1, vec2)
    return cp > 0 ? true : (cp == 0 ? dot2d(vec1, vec2) > 0 : false)
end
aligned_right_turn(o, a, b) = aligned_right_turn(triplet_to_vectors(o, a, b)...)
function aligned_right_turn(vec1, vec2)
    cp = cross2d(vec1, vec2)
    return cp < 0 ? true : (cp == 0 ? dot2d(vec1, vec2) > 0 : false)
end

# For the methods below, colinear vectors will only yield `true` if they are aligned (dot product < 0)
misaligned_left_turn(o, a, b) = misaligned_left_turn(triplet_to_vectors(o, a, b)...)
function misaligned_left_turn(vec1, vec2)
    cp = cross2d(vec1, vec2)
    return cp > 0 ? true : (cp == 0 ? 0 > dot2d(vec1, vec2) : false)
end
misaligned_right_turn(o, a, b) = misaligned_right_turn(triplet_to_vectors(o, a, b)...)
function misaligned_right_turn(vec1, vec2)
    cp = cross2d(vec1, vec2)
    return cp < 0 ? true : (cp == 0 ? 0 > dot2d(vec1, vec2) : false)
end