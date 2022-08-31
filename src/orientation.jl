@enum HullOrientation CW CCW
Base.:(!)(o::HullOrientation) = o === CCW ? CW : CCW

# default previous edge direction initializations
const UP = (0,1)
const DOWN = (0,-1)

# The following three methods are only meant to be called if the cross product is zero. They check if the next candidate edge is aligned with the previous edge.
isaligned(prevedge,nextedge) = dot2d(prevedge, nextedge) > 0
function isalignedfurther(orientation,prevedge, oa, ob) 
    dot_oa_ob = dot2d(oa, ob)
    if dot_oa_ob == 0           # if one of the vectors is of zero length...
        return sum(abs2, oa) == 0     # return true if OA is zero
    elseif dot_oa_ob > 0                    # if the two candidate vectors are pointing in the same direction...
        return dot_oa_ob > dot2d(oa, oa)    # return true if OB is longer
    else
        cp_prev_ob = cross2d(prevedge, ob)
        if cp_prev_ob == 0                  # if the previous edge, OA, and OB are all collinear, but OA and OB are not aligned,
            return isaligned(prevedge, ob)  # return true OB if it is aligned with the previous edge
        else                                                    # otherwise 
            return isorientedturn(orientation, prevedge, ob)    # return true if OB results in the proper orientation
        end
    end
end
function isalignedcloser(orientation,prevedge, oa, ob) 
    dot_oa_ob = dot2d(oa, ob)
    if dot_oa_ob == 0           # if one of the vectors is of zero length...
        return sum(abs2, oa) == 0     # return true if OA is zero
    elseif dot_oa_ob > 0                    # if the two candidate vectors are pointing in the same direction...
        return dot2d(oa, oa) > dot_oa_ob    # return true if OB is shorter
    else
        cp_prev_ob = cross2d(prevedge, ob)
        if cp_prev_ob == 0                  # if the previous edge, OA, and OB are all collinear, but OA and OB are not aligned,
            return isaligned(prevedge, ob)  # return true OB if it is aligned with the previous edge
        else                                                    # otherwise 
            return isorientedturn(orientation, prevedge, ob)    # return true if OB results in the proper orientation
        end
    end
end

# determine if the change from one vector to another represents a valid turn with the desired orientation and condition
function isvalidturn(orientation::HullOrientation, condition::Function, oa, ob)
    cp = cross2d(oa, ob)
    oriented = orientation === CCW ? (cp >= 0) : (cp <= 0)
    return oriented && condition(cp, oa, ob)
end
isvalidturn(orientation::HullOrientation, condition::Function, o, a, b) = isvalidturn(orientation, condition, (a[1] - o[1], a[2] - o[2]), (b[1] - o[1], b[2] - o[2]))

# colinear vectors will always yield `true`
isorientedturn(orientation::HullOrientation, oa, ob) = isvalidturn(orientation, (cp,v1,v2)->true, oa, ob)
isorientedturn(orientation::HullOrientation, o, a, b) = isorientedturn(orientation, (a[1] - o[1], a[2] - o[2]), (b[1] - o[1], b[2] - o[2]))


# colinear vectors will only yield `true` if they are aligned (dot product > 0). This is typically only useful with sorted points
isalignedturn(orientation, oa, ob) = isvalidturn(orientation, (cp, oa, ob) -> (cp != 0 ? true : isaligned(oa,ob)), oa, ob)
isalignedturn(orientation, o, a, b) = isalignedturn(orientation, (a[1] - o[1], a[2] - o[2]), (b[1] - o[1], b[2] - o[2]))


# colinear vectors will only yield `true` if the second is shorter than the first. This is typically only useful with sorted points
isshorterturn(orientation, oa, ob) = isvalidturn(orientation, (cp, oa, ob) -> (cp != 0 ? true : sum(abs2,oa) > sum(abs2,ob)), oa, ob)
isshorterturn(orientation, o, a, b) = isshorterturn(orientation, (a[1] - o[1], a[2] - o[2]), (b[1] - o[1], b[2] - o[2]))


# colinear vectors will only yield 'true' if they are aligned and if the second vector moves less than the first (relative to the direction of the previous edge)
iscloserturn(orientation, prevedge, oa, ob) = isvalidturn(orientation, (cp, oa, ob) -> (cp != 0 ? true : isalignedcloser(orientation,prevedge,oa,ob)), oa, ob)
iscloserturn(orientation, prevedge, o, a, b) = iscloserturn(orientation, prevedge, (a[1] - o[1], a[2] - o[2]), (b[1] - o[1], b[2] - o[2]))

# colinear vectors will only yield 'true' if they are aligned and if the second vector moves further than the first (relative to the direction of the previous edge)
isfurtherturn(orientation, prevedge, oa, ob) = isvalidturn(orientation, (cp, oa, ob) -> (cp != 0 ? true : isalignedfurther(orientation,prevedge,oa,ob)), oa, ob)
isfurtherturn(orientation, prevedge, o, a, b) = isfurtherturn(orientation, prevedge, (a[1] - o[1], a[2] - o[2]), (b[1] - o[1], b[2] - o[2]))
