@testset "generic axes" begin
    using OffsetArrays

    pts = [(i, j) for i in 1:5 for j in 1:5]

    # shifted-axes wrapper (sharpest test: exercises non-1-based indexing)
    opts = OffsetArray(pts, -10)
    # lazy wrapper (catches Array-ness / contiguity assumptions)
    vpts = view(pts, :)

    for (label, input) in (("offset", opts), ("view", vpts))
        @testset "$label input" begin
            @test monotonechain(input)       == monotonechain(pts)
            @test lower_monotonechain(input) == lower_monotonechain(pts)
            @test upper_monotonechain(input) == upper_monotonechain(pts)
            @test jarvismarch(input)         == jarvismarch(pts)
            @test lower_jarvismarch(input)   == lower_jarvismarch(pts)
            @test upper_jarvismarch(input)   == upper_jarvismarch(pts)

            let ref = MutableConvexHull{eltype(pts)}()
                mergepoints!(ref, pts)
                h = MutableConvexHull{eltype(pts)}()
                mergepoints!(h, input)
                @test h == ref
            end
            let ref = MutableLowerConvexHull{eltype(pts)}()
                mergepoints!(ref, pts)
                h = MutableLowerConvexHull{eltype(pts)}()
                mergepoints!(h, input)
                @test h == ref
            end
            let ref = MutableUpperConvexHull{eltype(pts)}()
                mergepoints!(ref, pts)
                h = MutableUpperConvexHull{eltype(pts)}()
                mergepoints!(h, input)
                @test h == ref
            end
        end
    end
end
