@testset "convex hulls, unique points" begin
    boxcoords = [(i,j) for i in 1:10 for j in 1:10]
    by = x -> (x[1], -x[2])
    n = 10
    
    hulltestset("lower convex hull", n, by, boxcoords, MutableLowerConvexHull, lower_jarvismarch)
    hulltestset("upper convex hull", n, by, boxcoords, MutableUpperConvexHull, upper_jarvismarch)
    hulltestset("convex hull",       n, by, boxcoords, MutableConvexHull,      jarvismarch)
end

@testset "convex hulls, duplicate points" begin
    boxcoords = [(i,j) for i in 1:10 for j in 1:10]
    boxcoords = [boxcoords..., boxcoords..., boxcoords...]
    by = x -> (x[1], -x[2])
    n = 10
    
    hulltestset("lower convex hull", n, by, boxcoords, MutableLowerConvexHull, lower_jarvismarch)
    hulltestset("upper convex hull", n, by, boxcoords, MutableUpperConvexHull, upper_jarvismarch)
    hulltestset("convex hull",       n, by, boxcoords, MutableConvexHull,      jarvismarch)
end

@testset "convex hulls, random data with duplicates" begin
    coords = [(randn(),randn()) for i in 1:10 for j in 1:10]
    by = x -> (x[1], -x[2])
    n = 10
    
    hulltestset("lower convex hull", n, by, coords, MutableLowerConvexHull, lower_jarvismarch)
    hulltestset("upper convex hull", n, by, coords, MutableUpperConvexHull, upper_jarvismarch)
    hulltestset("convex hull",       n, by, coords, MutableConvexHull,      jarvismarch)
end
@testset "hash, copy, and mergehulls" begin
    coords = [(randn(), randn()) for _ in 1:50]
    for H in (MutableConvexHull, MutableLowerConvexHull, MutableUpperConvexHull)
        @testset "$H" begin
            h = H{eltype(coords)}()
            for p in coords
                addpoint!(h, p)
            end

            # hash is consistent with `==`: equal hulls hash equally
            h2 = H{eltype(coords)}()
            for p in coords
                addpoint!(h2, p)
            end
            @test h == h2
            @test hash(h) == hash(h2)
            # usable as Set elements / Dict keys
            @test length(Set([h, h2])) == 1
            @test (Dict(h => 1)[h2]) == 1

            # copy is equal and independent
            hc = copy(h)
            @test hc == h
            @test hash(hc) == hash(h)
            @test hc isa H
            # no linked-list nodes are shared between original and copy
            @test !any(n1 === n2 for n1 in MCH.PointNodeIterator(h) for n2 in MCH.PointNodeIterator(hc))
            # mutating the copy leaves the original untouched, and vice versa
            origverts = collect(h)
            removepoint!(hc, hc.points.head.next)
            @test collect(h) == origverts
            copyverts = collect(hc)
            removepoint!(h, h.points.head.next)
            @test collect(hc) == copyverts

            # copying an empty hull yields an empty, independent hull
            empt = H{eltype(coords)}()
            ec = copy(empt)
            @test isempty(ec)
            @test ec == empt
        end
    end

    # mergehulls returns the merged hull without mutating its arguments
    @testset "mergehulls is non-mutating" begin
        a = MutableConvexHull{Tuple{Float64,Float64}}()
        b = MutableConvexHull{Tuple{Float64,Float64}}()
        for p in [(0.0,0.0),(1.0,0.0),(1.0,1.0),(0.0,1.0)]
            addpoint!(a, p)
        end
        for p in [(2.0,2.0),(3.0,2.0),(3.0,3.0),(2.0,3.0)]
            addpoint!(b, p)
        end
        averts = collect(a)
        bverts = collect(b)
        # `mergehulls` is exported, so it is callable unqualified
        @test :mergehulls in names(MutableConvexHulls)
        m = mergehulls(a, b)
        @test collect(a) == averts
        @test collect(b) == bverts
        merged = mergehulls!(copy(a), copy(b))
        @test m == merged
    end
end

@testset "removepoint! by value" begin
    boxcoords = [(i, j) for i in 1:10 for j in 1:10]
    coords = [boxcoords..., boxcoords...]   # include duplicate points
    for (H, truthfun) in ((MutableLowerConvexHull, lower_jarvismarch),
                          (MutableUpperConvexHull, upper_jarvismarch),
                          (MutableConvexHull, jarvismarch))
        @testset "$H" begin
            remaining = shuffle(coords)
            h = H{eltype(coords)}()
            mergepoints!(h, copy(remaining))
            for _ in 1:length(coords)
                p = rand(remaining)
                # drop one instance from the reference multiset, then by value from the hull
                deleteat!(remaining, findfirst(==(p), remaining))
                removepoint!(h, p)
                isempty(remaining) ? (@test isempty(h)) : (@test h == truthfun(remaining))
            end
            # fail-fast when no contained point equals the value
            @test_throws ArgumentError removepoint!(h, (-1, -1))
        end
    end
end

@testset "node iterator eltype" begin
    coords = [(i, j) for i in 1:5 for j in 1:5]
    for H in (MutableLowerConvexHull, MutableUpperConvexHull, MutableConvexHull)
        @testset "$H" begin
            h = H{eltype(coords)}()
            mergepoints!(h, coords)

            hulliter = MCH.HullNodeIterator(h)
            pointiter = MCH.PointNodeIterator(h)

            # eltype reports the concrete node type rather than `Any`
            @test eltype(hulliter) === MCH.nodetype(h.hull) !== Any
            @test eltype(pointiter) === MCH.nodetype(h.points) !== Any
            @test Base.IteratorEltype(hulliter) === Base.HasEltype()
            @test Base.IteratorEltype(pointiter) === Base.HasEltype()

            # collecting yields a vector typed by that eltype
            @test collect(hulliter) isa Vector{eltype(hulliter)}
            @test collect(pointiter) isa Vector{eltype(pointiter)}
        end
    end
end
