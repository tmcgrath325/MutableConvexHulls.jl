@testset "chan hulls, unique points" begin
    boxcoords = [(i,j) for i in 1:10 for j in 1:10]
    by = x -> (x[1], -x[2])
    n = 10
    
    chanhulltestset("lower chan hull", n, by, boxcoords, ChanLowerConvexHull, lower_jarvismarch)
    chanhulltestset("upper chan hull", n, by, boxcoords, ChanUpperConvexHull, upper_jarvismarch)
    chanhulltestset("chan hull",       n, by, boxcoords, ChanConvexHull,      jarvismarch)
end

@testset "chan hulls, duplicate points" begin
    boxcoords = [(i,j) for i in 1:10 for j in 1:10]
    boxcoords = [boxcoords..., boxcoords..., boxcoords...]
    by = x -> (x[1], -x[2])
    n = 10
    
    chanhulltestset("lower chan hull", n, by, boxcoords, ChanLowerConvexHull, lower_jarvismarch)
    chanhulltestset("upper chan hull", n, by, boxcoords, ChanUpperConvexHull, upper_jarvismarch)
    chanhulltestset("chan hull",       n, by, boxcoords, ChanConvexHull,      jarvismarch)
end

@testset "chan hulls, random data with duplicates" begin
    coords = [(randn(),randn()) for i in 1:10 for j in 1:10]
    coords = [coords..., coords..., coords...]
    by = x -> (x[1], -x[2])
    n = 10

    chanhulltestset("lower chan hull", n, by, coords, ChanLowerConvexHull, lower_jarvismarch)
    chanhulltestset("upper chan hull", n, by, coords, ChanUpperConvexHull, upper_jarvismarch)
    chanhulltestset("chan hull",       n, by, coords, ChanConvexHull,      jarvismarch)
end

@testset "chan hull AbstractConvexHull interface" begin
    coords = [(i,j) for i in 1:10 for j in 1:10]
    for H in (ChanConvexHull, ChanLowerConvexHull, ChanUpperConvexHull)
        @testset "$H" begin
            h = H{eltype(coords)}()
            mergepoints!(h, coords)
            # Chan hulls participate in the AbstractConvexHull interface
            @test h isa AbstractConvexHull{eltype(coords)}
            @test length(h) > 0
            @test eltype(h) === eltype(coords)
            @test !isempty(h)
            # show iterates the hull vertices; it must not touch a `points` field (Chan hulls lack one)
            @test occursin(string(typeof(h)), sprint(show, h))
            # an independently built hull over the same points compares equal
            h2 = H{eltype(coords)}()
            mergepoints!(h2, coords)
            @test h == h2
            # empty(h) yields an empty hull preserving attributes; h is untouched
            he = empty(h)
            @test he isa H
            @test isempty(he)
            @test length(h) > 0
            @test he.orientation == h.orientation && he.collinear == h.collinear && he.sortedby == h.sortedby
            # empty!(h) clears the hull in place and resets to a single subhull
            empty!(h)
            @test isempty(h)
            @test length(h.subhulls) == 1
            @test isempty(h.subhulls[1])
        end
    end
end