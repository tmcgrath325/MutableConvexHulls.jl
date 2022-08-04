using MutableConvexHulls: oriented_turn, aligned_turn, closer_turn, further_turn

@testset "orientation" begin
    # set up a simple grid of coordinates for testing
    boxcoords = [(i,j) for i in 1:4 for j in 0:3]
    # use angles in the x-y plane to confirm orientation
    wrapangle(angle) = angle - sign(angle) * 2π * ceil((abs(angle) - π)/(2π)) # always [-π,π]

    @testset "basic left/right" begin
        for o in boxcoords
            for a in boxcoords
                for b in boxcoords
                    left_result = oriented_turn(CCW, o, a, b)
                    right_result = oriented_turn(CW, o, a, b)
                    oa = a .- o
                    ob = b .- o
                    @test left_result == oriented_turn(CCW, oa, ob)
                    @test right_result == oriented_turn(CW, oa, ob)
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
                    left_result = aligned_turn(CCW, o, a, b)
                    right_result = aligned_turn(CW, o, a, b)
                    oa = a .- o
                    ob = b .- o
                    @test left_result == aligned_turn(CCW, oa, ob)
                    @test right_result == aligned_turn(CW, oa, ob)
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
        for o in boxcoords
            for a in boxcoords
                for b in boxcoords
                    closer_left_result = closer_turn(CCW, o, a, b)
                    closer_right_result = closer_turn(CW, o, a, b)
                    further_left_result = further_turn(CCW, o, a, b)
                    further_right_result = further_turn(CW, o, a, b)
                    oa = a .- o
                    ob = b .- o
                    @test closer_left_result == closer_turn(CCW, oa, ob)
                    @test closer_right_result == closer_turn(CW, oa, ob)
                    @test further_left_result == further_turn(CCW, oa, ob)
                    @test further_right_result == further_turn(CW, oa, ob)
                    angle1 = atan(oa[2], oa[1])
                    angle2 = atan(ob[2], ob[1])
                    anglediff = wrapangle(angle2 - angle1)
                    if sum(abs2, oa) == 0 || sum(abs2, ob) == 0 || abs(anglediff) ≈ π
                        @test !closer_left_result && !closer_right_result && !further_left_result && !further_right_result
                    elseif abs(anglediff) ≈ 0
                        lensq1 = sum(abs2, oa)
                        lensq2 = sum(abs2, ob)
                        @test lensq1 > lensq2 ? (closer_left_result && closer_right_result) : (!closer_left_result && !closer_right_result)
                        @test lensq2 > lensq1 ? (further_left_result && further_right_result) : (!further_left_result && !further_right_result)
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