@enum HullOrientation CW CCW

"Clockwise hull orientation: hull vertices are listed in clockwise order. See also [`CCW`](@ref)."
CW

"Counterclockwise hull orientation: hull vertices are listed in counterclockwise order. See also [`CW`](@ref)."
CCW

Base.:(!)(o::HullOrientation) = o === CCW ? CW : CCW

# default previous edge direction initializations
const UP = (0, 1)
const DOWN = (0, -1)

# The following three methods are only meant to be called if the cross product is zero. They check if the next candidate edge is aligned with the previous edge.
isaligned(prevedge, nextedge) = dot2d(prevedge, nextedge) > 0
function isalignedfurther(orientation, prevedge, oa, ob, ab)
    dot_oa_ob = dot2d(oa, ob)
    if dot_oa_ob == 0           # if one of the vectors is of zero length...
        return sum(abs2, oa) == 0     # return true if OA is zero
    elseif dot_oa_ob > 0                    # if the two candidate vectors are pointing in the same direction...
        return dot_oa_ob >= dot2d(oa, oa)   # return true if OB is longer (or equal)
    else
        cp_prev_ob = cross2d(prevedge, ob)
        if cp_prev_ob == 0                  # if the previous edge, OA, and OB are all collinear, but OA and OB are not aligned,
            return isaligned(prevedge, ob)  # return true OB if it is aligned with the previous edge
        else                                                    # otherwise
            return isorientedturn_vec(orientation, prevedge, ob, ab)    # return true if OB results in the proper orientation
        end
    end
end
function isalignedcloser(orientation, prevedge, oa, ob, ab)
    dot_oa_ob = dot2d(oa, ob)
    if dot_oa_ob == 0           # if one of the vectors is of zero length...
        return sum(abs2, oa) == 0     # return true if OA is zero
    elseif dot_oa_ob > 0                    # if the two candidate vectors are pointing in the same direction...
        return dot2d(oa, oa) >= dot_oa_ob   # return true if OB is shorter (or equal)
    else
        cp_prev_ob = cross2d(prevedge, ob)
        if cp_prev_ob == 0                  # if the previous edge, OA, and OB are all collinear, but OA and OB are not aligned,
            return isaligned(prevedge, ob)  # return true OB if it is aligned with the previous edge
        else                                                    # otherwise
            return isorientedturn_vec(orientation, prevedge, ob, ab)    # return true if OB results in the proper orientation
        end
    end
end

# determine if the change from one vector to another represents a valid turn with the desired orientation and condition
function isvalidturn(orientation::HullOrientation, condition::Function, oa, ob, ab)
    cp = cross2d(oa, ab)
    oriented = orientation === CCW ? (cp >= 0) : (cp <= 0)
    return oriented && condition(cp, oa, ob, ab)
end

# Float64 fast path for the point-based turn predicates below. `fastturn`
# evaluates cross2d(a - o, b - a) in Float64 and certifies its sign against a
# bound on the worst-case rounding error: with unit roundoff u = eps()/2, each
# coordinate difference carries relative error at most u and each product at
# most 3u, so the absolute error stays below 4u * (|l| + |r|) up to O(u^2)
# terms; the factor 8 * eps() = 16u leaves a 4x margin. A certified nonzero
# sign decides every turn predicate outright, because their collinearity
# conditions only apply when the cross product is exactly zero. Uncertain cases
# (ties and near-ties) return `nothing` and must be recomputed with
# DoubleFloats via `sub2d`.
const CROSSERRFACTOR = 8 * eps(Float64)
# The relative-error analysis requires both products to be in the normal
# range; below this floor one of them may be subnormal, voiding the bound, so
# certification is refused.
const CROSSNORMALFLOOR = 0x1p-1020

# Float32 converts to Float64 exactly, so certification remains valid for the
# original coordinates.
function fastturn(orientation::HullOrientation, o::P, a::P, b::P) where {P <: Union{NTuple{2, Float32}, NTuple{2, Float64}}}
    oax = Float64(a[1]) - Float64(o[1])
    oay = Float64(a[2]) - Float64(o[2])
    abx = Float64(b[1]) - Float64(a[1])
    aby = Float64(b[2]) - Float64(a[2])
    l = oax * aby
    r = oay * abx
    det = l - r
    magnitude = abs(l) + abs(r)
    magnitude < CROSSNORMALFLOOR && return nothing
    err = CROSSERRFACTOR * magnitude
    det > err && return orientation === CCW
    det < -err && return orientation === CW
    return nothing  # includes non-finite inputs, for which both comparisons are false
end
# Other coordinate types (integers can overflow, exotic reals have their own
# rounding behavior) always take the DoubleFloat path.
fastturn(orientation::HullOrientation, o, a, b) = nothing

# colinear vectors will always yield `true`
isorientedturn_vec(orientation::HullOrientation, oa, ob, ab) = isvalidturn(orientation, (cp, v1, v2, v3) -> true, oa, ob, ab)
function isorientedturn(orientation::HullOrientation, o, a, b)
    f = fastturn(orientation, o, a, b)
    f === nothing || return f
    return isorientedturn_vec(orientation, sub2d(a, o), sub2d(b, o), sub2d(b, a))
end

# colinear vectors will only yield `true` if they are aligned (dot product > 0). This is typically only useful with sorted points
isalignedturn_vec(orientation, oa, ob, ab) = isvalidturn(orientation, (cp, oa, ob, ab) -> (cp != 0 ? true : isaligned(oa, ob)), oa, ob, ab)
function isalignedturn(orientation, o, a, b)
    f = fastturn(orientation, o, a, b)
    f === nothing || return f
    return isalignedturn_vec(orientation, sub2d(a, o), sub2d(b, o), sub2d(b, a))
end

# colinear vectors will only yield `true` if the second is shorter than the first. This is typically only useful with sorted points
isshorterturn_vec(orientation, oa, ob, ab) = isvalidturn(orientation, (cp, oa, ob, ab) -> (cp != 0 ? true : sum(abs2, oa) > sum(abs2, ob)), oa, ob, ab)
function isshorterturn(orientation, o, a, b)
    f = fastturn(orientation, o, a, b)
    f === nothing || return f
    return isshorterturn_vec(orientation, sub2d(a, o), sub2d(b, o), sub2d(b, a))
end

# colinear vectors will only yield 'true' if they are aligned and if the second vector moves less than the first (relative to the direction of the previous edge)
iscloserturn_vec(orientation, prevedge, oa, ob, ab) = isvalidturn(orientation, (cp, oa, ob, ab) -> (cp != 0 ? true : isalignedcloser(orientation, prevedge, oa, ob, ab)), oa, ob, ab)
function iscloserturn(orientation, prevedge, o, a, b)
    f = fastturn(orientation, o, a, b)
    f === nothing || return f
    return iscloserturn_vec(orientation, prevedge, sub2d(a, o), sub2d(b, o), sub2d(b, a))
end

# colinear vectors will only yield 'true' if they are aligned and if the second vector moves further than the first (relative to the direction of the previous edge)
isfurtherturn_vec(orientation, prevedge, oa, ob, ab) = isvalidturn(orientation, (cp, oa, ob, ab) -> (cp != 0 ? true : isalignedfurther(orientation, prevedge, oa, ob, ab)), oa, ob, ab)
function isfurtherturn(orientation, prevedge, o, a, b)
    f = fastturn(orientation, o, a, b)
    f === nothing || return f
    return isfurtherturn_vec(orientation, prevedge, sub2d(a, o), sub2d(b, o), sub2d(b, a))
end
