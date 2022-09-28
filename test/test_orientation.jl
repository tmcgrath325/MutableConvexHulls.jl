using MutableConvexHulls: isaligned, cross2d
using MutableConvexHulls: isorientedturn, isalignedturn, iscloserturn, isfurtherturn
using MutableConvexHulls: isorientedturn_vec, isalignedturn_vec, iscloserturn_vec, isfurtherturn_vec

@testset "orientation" begin
    # set up a simple grid of coordinates for testing
    boxcoords = [(i,j) for i in 1:4 for j in 0:3]
    # use angles in the x-y plane to confirm orientation
    wrapangle(angle) = angle - sign(angle) * 2π * ceil((abs(angle) - π)/(2π)) # always [-π,π]

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
        prevdirections = filter(x->x!=(0,0), [(i,j) for i in (-1,0,1) for j in (-1,0,1)])
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
                            # @show prev, oa, ob
                            # sleep(0.1)
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
end