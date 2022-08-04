# clockwise/counterclockwise selection
struct CW end
struct CCW end
const HullOrientation = Type{<:Union{CW,CCW}}

Base.:(!)(::Type{CW}) = CCW
Base.:(!)(::Type{CCW}) = CW

aligned(oa, ob) = dot2d(oa, ob) > 0
aligned_further(oa, ob) = dot2d(oa, ob) > dot2d(oa, oa) > 0
aligned_closer(oa, ob) = dot2d(oa, oa) > dot2d(oa, ob) > 0

# determine if the change from one vector to another represents a valid turn with the desired orientation and condition
function valid_turn(orientation::HullOrientation, condition::Function, oa, ob)
    cp = cross2d(oa, ob)
    oriented = orientation === CCW ? (cp >= 0) : (cp <= 0)
    return oriented && condition(cp, oa, ob)
end
valid_turn(orientation::HullOrientation, condition::Function, o, a, b) = valid_turn(orientation, condition, a .- o, b .- o)

# colinear vectors will always yield `true`
oriented_turn(orientation::HullOrientation, args...) = valid_turn(orientation, (cp,v1,v2)->true, args...)

# colinear vectors will only yield `true` if they are aligned (dot product > 0)
aligned_turn(orientation, args...) = valid_turn(orientation, (cp, oa, ob) -> (cp != 0 ? true : aligned(oa,ob)), args...)

# colinear vectors will only yield 'true' if they are aligned and if the second vector is shorter than the first
closer_turn(orientation, args...) = valid_turn(orientation, (cp, oa, ob) -> (cp != 0 ? true : aligned_closer(oa,ob)), args...)

# colinear vectors will only yield 'true' if they are aligned and if the second vector is longer than the first
further_turn(orientation, args...) = valid_turn(orientation, (cp, oa, ob) -> (cp != 0 ? true : aligned_further(oa,ob)), args...)
