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
@testset "chan hull copy" begin
    coords = [(i + 0.5*j, j - 0.3*i) for i in 1:10 for j in 1:10]
    for H in (ChanConvexHull, ChanLowerConvexHull, ChanUpperConvexHull)
        @testset "$H" begin
            h = H{eltype(coords)}()
            mergepoints!(h, coords)

            hc = copy(h)
            @test hc isa H
            @test hc == h
            @test hash(hc) == hash(h)
            # the copy owns independent subhulls
            @test all(s1 !== s2 for s1 in h.subhulls for s2 in hc.subhulls)
            # mutating the copy leaves the original untouched
            origverts = collect(h)
            addpoint!(hc, (100.0, 100.0))
            @test collect(h) == origverts
        end
    end
end

@testset "chan removepoint! by HullNode" begin
    coords = [(i, j) for i in 1:5 for j in 1:5]
    for (H, truthfun) in ((ChanConvexHull,      jarvismarch),
                          (ChanLowerConvexHull,  lower_jarvismarch),
                          (ChanUpperConvexHull,  upper_jarvismarch))
        @testset "$H" begin
            h = H{eltype(coords)}()
            mergepoints!(h, coords)
            # HullNode comes from a subhull's hull list, not from h.hull
            subhull = h.subhulls[1]
            hullnode = subhull.hull.head.next
            vertex = hullnode.data

            h2 = H{eltype(coords)}()
            mergepoints!(h2, coords)

            removepoint!(h, hullnode)  # by HullNode (chanhull.jl dispatch)
            removepoint!(h2, vertex)   # by value — same expected result
            @test h == h2
        end
    end
end

@testset "chan removepoint! by value and insidehull" begin
    boxcoords = [(i, j) for i in 1:10 for j in 1:10]
    coords = [boxcoords..., boxcoords...]   # include duplicate points
    queries = [(i, j) for i in 0:11 for j in 0:11]   # interior, boundary, and outside
    for (H, R, truthfun) in ((ChanLowerConvexHull, MutableLowerConvexHull, lower_jarvismarch),
                             (ChanUpperConvexHull, MutableUpperConvexHull, upper_jarvismarch),
                             (ChanConvexHull,      MutableConvexHull,      jarvismarch))
        @testset "$H" begin
            # insidehull on a Chan hull agrees with the equivalent regular hull
            hc = H{eltype(coords)}(); mergepoints!(hc, copy(coords))
            hr = R{eltype(coords)}(); mergepoints!(hr, copy(coords))
            @test all(insidehull(q, hc) == insidehull(q, hr) for q in queries)

            # removepoint! by value tracks the truth hull of the remaining points
            remaining = shuffle(coords)
            h = H{eltype(coords)}()
            mergepoints!(h, copy(remaining))
            for _ in 1:length(coords)
                p = rand(remaining)
                deleteat!(remaining, findfirst(==(p), remaining))
                removepoint!(h, p)
                isempty(remaining) ? (@test isempty(h)) : (@test h == truthfun(remaining))
            end
            # fail-fast when no contained point equals the value
            @test_throws ArgumentError removepoint!(h, (-1, -1))
        end
    end
end
