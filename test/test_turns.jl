using MutableConvexHulls: left_turn, right_turn
using MutableConvexHulls: aligned_left_turn, aligned_right_turn
using MutableConvexHulls: misaligned_left_turn, misaligned_right_turn

@testset "turns" begin
    boxcoords = [(i,j) for i in 1:4 for j in 0:3]
    wrapangle(angle) = angle - sign(angle) * 2π * ceil((abs(angle) - π)/(2π)) # always [-π,π]

    @testset "basic" begin
        for c1 in boxcoords
            for c2 in boxcoords
                for c3 in boxcoords
                    leftresult = left_turn(c1, c2, c3)
                    rightresult = right_turn(c1, c2, c3)
                    vec1 = c2 .- c1
                    vec2 = c3 .- c1
                    @test leftresult == left_turn(vec1, vec2)
                    @test rightresult == right_turn(vec1, vec2)
                    angle1 = atan(vec1[2], vec1[1])
                    angle2 = atan(vec2[2], vec2[1])
                    anglediff = wrapangle(angle2 - angle1)
                    if sum(abs2, vec1) == 0 || sum(abs2, vec2) == 0 || anglediff == 0 || abs(anglediff) ≈ π     # angles are in increments of π/4 so we don't need to worry about floating point error
                        @test leftresult && rightresult
                    else 
                        @test leftresult == (wrapangle(anglediff) >= 0)     # if a left turn, the angle difference should be [0, π]
                        @test rightresult == (wrapangle(anglediff) <= 0)    # if a right turn, the angle difference should be [-π, 0]
                        @test leftresult != rightresult
                    end
                end
            end
        end
    end

    @testset "aligned" begin
        for c1 in boxcoords
            for c2 in boxcoords
                for c3 in boxcoords
                    leftresult = aligned_left_turn(c1, c2, c3)
                    rightresult = aligned_right_turn(c1, c2, c3)
                    vec1 = c2 .- c1
                    vec2 = c3 .- c1
                    @test leftresult == aligned_left_turn(vec1, vec2)
                    @test rightresult == aligned_right_turn(vec1, vec2)
                    angle1 = atan(vec1[2], vec1[1])
                    angle2 = atan(vec2[2], vec2[1])
                    anglediff = wrapangle(angle2 - angle1)
                    if sum(abs2, vec1) == 0 || sum(abs2, vec2) == 0 || abs(anglediff) ≈ π
                        @test !leftresult && !rightresult
                    elseif anglediff == 0
                        @test leftresult && rightresult
                    else 
                        @test leftresult == (wrapangle(anglediff) >= 0)
                        @test rightresult == (wrapangle(anglediff) <= 0)
                        @test leftresult != rightresult
                    end
                end
            end
        end
    end

    @testset "misaligned" begin
        for c1 in boxcoords
            for c2 in boxcoords
                for c3 in boxcoords
                    leftresult = misaligned_left_turn(c1, c2, c3)
                    rightresult = misaligned_right_turn(c1, c2, c3)
                    vec1 = c2 .- c1
                    vec2 = c3 .- c1
                    @test leftresult == misaligned_left_turn(vec1, vec2)
                    @test rightresult == misaligned_right_turn(vec1, vec2)
                    angle1 = atan(vec1[2], vec1[1])
                    angle2 = atan(vec2[2], vec2[1])
                    anglediff = wrapangle(angle2 - angle1)
                    if sum(abs2, vec1) == 0 || sum(abs2, vec2) == 0 || anglediff == 0
                        @test !leftresult && !rightresult
                    elseif abs(anglediff) ≈ π
                        @test leftresult && rightresult
                    else 
                        @test leftresult == (wrapangle(anglediff) >= 0)
                        @test rightresult == (wrapangle(anglediff) <= 0)
                        @test leftresult != rightresult
                    end
                end
            end
        end
    end
end