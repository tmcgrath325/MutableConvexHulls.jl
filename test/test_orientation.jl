using MutableConvexHulls: isaligned, cross2d, fastturn
using MutableConvexHulls: isorientedturn, isalignedturn, isshorterturn, iscloserturn, isfurtherturn
using MutableConvexHulls: isorientedturn_vec, isalignedturn_vec, iscloserturn_vec, isfurtherturn_vec

@testset "orientation" begin
    # set up a simple grid of coordinates for testing
    boxcoords = [(i, j) for i in 1:4 for j in 0:3]
    # use angles in the x-y plane to confirm orientation
    wrapangle(angle) = angle - sign(angle) * 2π * ceil((abs(angle) - π) / (2π)) # always [-π,π]

    @testset "basic left/right" begin
        for o in boxcoords
            for a in boxcoords
                for b in boxcoords
                    left_result = isorientedturn(CCW, o, a, b)
                    right_result = isorientedturn(CW, o, a, b)
                    oa = a .- o
                    ob = b .- o
                    ab = b .- a
                    @test left_result == isorientedturn_vec(CCW, oa, ob, ab)
                    @test right_result == isorientedturn_vec(CW, oa, ob, ab)
                    angle1 = atan(oa[2], oa[1])
                    angle2 = atan(ob[2], ob[1])
                    anglediff = wrapangle(angle2 - angle1)
                    if sum(abs2, oa) == 0 || sum(abs2, ob) == 0 || anglediff == 0 || abs(anglediff) ≈ π     # angles are in increments of π/4 so we don't need to worry about floating point error
                        @test left_result && right_result
                    else
                        @test left_result == (wrapangle(anglediff) > 0)     # if a left turn, the angle difference should be [0, π]
                        @test right_result == (wrapangle(anglediff) < 0)    # if a right turn, the angle difference should be [-π, 0]
                        @test left_result != right_result
                    end
                end
            end
        end
    end

    @testset "aligned" begin
        for o in boxcoords
            for a in boxcoords
                for b in boxcoords
                    left_result = isalignedturn(CCW, o, a, b)
                    right_result = isalignedturn(CW, o, a, b)
                    oa = a .- o
                    ob = b .- o
                    ab = b .- a
                    @test left_result == isalignedturn_vec(CCW, oa, ob, ab)
                    @test right_result == isalignedturn_vec(CW, oa, ob, ab)
                    angle1 = atan(oa[2], oa[1])
                    angle2 = atan(ob[2], ob[1])
                    anglediff = wrapangle(angle2 - angle1)
                    if sum(abs2, oa) == 0 || sum(abs2, ob) == 0 || abs(anglediff) ≈ π
                        @test !left_result && !right_result
                    elseif anglediff == 0
                        @test left_result && right_result
                    else
                        @test left_result == (wrapangle(anglediff) > 0)
                        @test right_result == (wrapangle(anglediff) < 0)
                        @test left_result != right_result
                    end
                end
            end
        end
    end

    @testset "closer/further" begin
        prevdirections = filter(x -> x != (0, 0), [(i, j) for i in (-1, 0, 1) for j in (-1, 0, 1)])
        for o in boxcoords
            for a in boxcoords
                for b in boxcoords
                    for prev in prevdirections
                        closer_left_result = iscloserturn(CCW, prev, o, a, b)
                        closer_right_result = iscloserturn(CW, prev, o, a, b)
                        further_left_result = isfurtherturn(CCW, prev, o, a, b)
                        further_right_result = isfurtherturn(CW, prev, o, a, b)
                        oa = a .- o
                        ob = b .- o
                        ab = b .- a
                        @test closer_left_result == iscloserturn_vec(CCW, prev, oa, ob, ab)
                        @test closer_right_result == iscloserturn_vec(CW, prev, oa, ob, ab)
                        @test further_left_result == isfurtherturn_vec(CCW, prev, oa, ob, ab)
                        @test further_right_result == isfurtherturn_vec(CW, prev, oa, ob, ab)
                        angle1 = atan(oa[2], oa[1])
                        angle2 = atan(ob[2], ob[1])
                        anglediff = wrapangle(angle2 - angle1)
                        if sum(abs2, oa) == 0
                            @test closer_left_result && closer_right_result && further_left_result && further_right_result
                        elseif sum(abs2, ob) == 0
                            @test !closer_left_result && !closer_right_result && !further_left_result && !further_right_result
                        elseif abs(anglediff) ≈ π
                            if cross2d(prev, oa) == 0
                                @test closer_left_result == further_left_result == closer_right_result == further_right_result == !isaligned(prev, oa)
                            else
                                @test closer_left_result == further_left_result == isorientedturn_vec(CCW, prev, ob, ab)
                                @test closer_right_result == further_right_result == isorientedturn_vec(CW, prev, ob, ab)
                            end
                        elseif abs(anglediff) ≈ 0
                            lensq1 = sum(abs2, oa)
                            lensq2 = sum(abs2, ob)
                            @test lensq1 >= lensq2 ? (closer_left_result && closer_right_result) : (!closer_left_result && !closer_right_result)
                            @test lensq2 >= lensq1 ? (further_left_result && further_right_result) : (!further_left_result && !further_right_result)
                        else
                            @test closer_left_result == further_left_result == (wrapangle(anglediff) > 0)
                            @test closer_right_result == further_right_result == (wrapangle(anglediff) < 0)
                            @test closer_left_result != closer_right_result
                        end
                    end
                end
            end
        end
    end

    @testset "Float64 fast path" begin
        # exact reference: sign of cross2d(a - o, b - a) in rational arithmetic
        function exactsign(o, a, b)
            R = Rational{BigInt}
            oa = (R(a[1]) - R(o[1]), R(a[2]) - R(o[2]))
            ab = (R(b[1]) - R(a[1]), R(b[2]) - R(a[2]))
            return Int(sign(oa[1] * ab[2] - oa[2] * ab[1]))
        end

        @testset "dispatch and guards" begin
            # only Float32/Float64 coordinates take the fast path
            @test fastturn(CCW, (0.0, 0.0), (1.0, 0.0), (1.0, 1.0)) === true
            @test fastturn(CW, (0.0, 0.0), (1.0, 0.0), (1.0, 1.0)) === false
            @test fastturn(CCW, (0.0f0, 0.0f0), (1.0f0, 0.0f0), (1.0f0, 1.0f0)) === true
            @test fastturn(CCW, (0, 0), (1, 0), (1, 1)) === nothing
            # exact ties cannot be certified
            @test fastturn(CCW, (0.0, 0.0), (1.0, 1.0), (2.0, 2.0)) === nothing
            # coordinates whose products may be subnormal refuse certification...
            @test fastturn(CCW, (0.0, 0.0), (1.0e-160, 0.0), (1.0e-160, 1.0e-160)) === nothing
            # ...but the DoubleFloat path still resolves the turn
            @test isorientedturn(CCW, (0.0, 0.0), (1.0e-160, 0.0), (1.0e-160, 1.0e-160))
            @test !isorientedturn(CW, (0.0, 0.0), (1.0e-160, 0.0), (1.0e-160, 1.0e-160))
            # non-finite input falls back rather than certifying
            @test fastturn(CCW, (0.0, 0.0), (Inf, 0.0), (1.0, 1.0)) === nothing
        end

        @testset "magnitude-disparate points" begin
            # a and b lie exactly on a line through the origin; o is a tiny offset
            # from that line which plain Float64 subtraction absorbs entirely
            # (2^-30 is below half an ulp of 2^40), so a naive Float64 cross
            # product is exactly zero. The predicates must resolve the strict turn.
            a = (2.0^40, 2.0^40)
            b = (2.0^41, 2.0^41)
            o = (2.0^-30, -(2.0^-30))       # exact cross2d(a - o, b - a) == -2^11
            @test fastturn(CW, o, a, b) === nothing
            @test isorientedturn(CW, o, a, b)
            @test !isorientedturn(CCW, o, a, b)
            omirror = (-(2.0^-30), 2.0^-30) # exact cross2d(a - o, b - a) == +2^11
            @test isorientedturn(CCW, omirror, a, b)
            @test !isorientedturn(CW, omirror, a, b)
        end

        @testset "agreement with exact arithmetic" begin
            rng = Random.Xoshiro(0x2026)
            triples = NTuple{3, NTuple{2, Float64}}[]
            for _ in 1:300  # huge nearly-collinear pair with a small off-line point
                a = (1.0e12 * rand(rng), 1.0e12 * rand(rng))
                push!(triples, ((1.0e-4 * randn(rng), 1.0e-4 * randn(rng)), a, 2.0 .* a))
            end
            for _ in 1:300  # borderline collinear at unit scale
                o = (rand(rng), rand(rng))
                d = (rand(rng) - 0.5, rand(rng) - 0.5)
                t, s = rand(rng), 1.0 + rand(rng)
                push!(triples, (o, o .+ t .* d, o .+ s .* d))
            end
            randmagnitude(rng) = rand(rng, (-1, 1)) * 10.0^(12 * rand(rng))
            for _ in 1:300  # coordinates spanning 12 orders of magnitude
                push!(triples, ((randmagnitude(rng), randmagnitude(rng)), (randmagnitude(rng), randmagnitude(rng)), (randmagnitude(rng), randmagnitude(rng))))
            end
            for _ in 1:300  # uniform
                push!(triples, ((rand(rng), rand(rng)), (rand(rng), rand(rng)), (rand(rng), rand(rng))))
            end
            prevedge = (1.0, 0.0)
            for (o, a, b) in triples
                s = exactsign(o, a, b)
                # certified fast-path signs must match exact arithmetic
                f = fastturn(CCW, o, a, b)
                if f !== nothing
                    @test f == (s > 0)
                end
                g = fastturn(CW, o, a, b)
                if g !== nothing
                    @test g == (s < 0)
                end
                @test isorientedturn(CCW, o, a, b) == (s >= 0)
                @test isorientedturn(CW, o, a, b) == (s <= 0)
                if s != 0   # with a strictly nonzero cross product, every turn predicate reduces to the sign test
                    for orient in (CCW, CW)
                        expected = orient === CCW ? s > 0 : s < 0
                        @test isalignedturn(orient, o, a, b) == expected
                        @test isshorterturn(orient, o, a, b) == expected
                        @test iscloserturn(orient, prevedge, o, a, b) == expected
                        @test isfurtherturn(orient, prevedge, o, a, b) == expected
                    end
                end
            end
        end
    end
end
