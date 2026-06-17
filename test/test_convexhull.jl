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

@testset "constructor keyword API" begin
    # All configuration arguments are keyword-only; positional calls now throw.
    T = Tuple{Float64,Float64}
    by = x -> x[1]
    for H in (MutableConvexHull, MutableLowerConvexHull, MutableUpperConvexHull)
        # Defaults
        h = H{T}()
        @test h.orientation === CCW && h.collinear === false && h.sortedby === identity
        # All keywords explicit
        h2 = H{T}(; orientation=CW, collinear=true, sortedby=by)
        @test h2.orientation === CW && h2.collinear === true && h2.sortedby === by
        # Positional call now throws
        @test_throws MethodError H{T}(CCW, false, identity)
    end
    for H in (ChanConvexHull, ChanLowerConvexHull, ChanUpperConvexHull)
        h = H{T}()
        @test h.orientation === CCW && h.collinear === false && h.sortedby === identity
        h2 = H{T}(; orientation=CW, collinear=true, sortedby=by)
        @test h2.orientation === CW && h2.collinear === true && h2.sortedby === by
        @test_throws MethodError H{T}(CCW, false, identity)
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

            # node-based HullNodeIterator constructor: starting from a specific node
            # produces the same suffix as collecting from that position in the full hull
            allnodes = collect(MCH.HullNodeIterator(h))
            @test collect(MCH.HullNodeIterator(allnodes[1])) == allnodes
            length(allnodes) >= 2 &&
                @test collect(MCH.HullNodeIterator(allnodes[2])) == allnodes[2:end]
        end
    end
end

@testset "removepoint! by HullNode" begin
    coords = [(i, j) for i in 1:5 for j in 1:5]
    for (H, truthfun) in ((MutableLowerConvexHull, lower_jarvismarch),
                          (MutableUpperConvexHull, upper_jarvismarch),
                          (MutableConvexHull,      jarvismarch))
        @testset "$H" begin
            h = H{eltype(coords)}()
            mergepoints!(h, coords)
            hullnode = h.hull.head.next
            vertex = hullnode.data

            h2 = H{eltype(coords)}()
            mergepoints!(h2, coords)

            removepoint!(h, hullnode)  # by HullNode
            removepoint!(h2, vertex)   # by value — same expected result
            @test h == h2
        end
    end
end

@testset "mergepoints! Matrix input" begin
    coords = [(i, j) for i in 1:5 for j in 1:5]
    m = [p[k] for p in coords, k in 1:2]
    for H in (MutableLowerConvexHull, MutableUpperConvexHull, MutableConvexHull)
        @testset "$H" begin
            h_vec = H{eltype(coords)}(); mergepoints!(h_vec, coords)
            h_mat = H{eltype(coords)}(); mergepoints!(h_mat, m)
            @test h_vec == h_mat
        end
    end
end

@testset "insidehull with AbstractNode" begin
    # The convex hull of the 5×5 grid is the corner rectangle.
    # Pass hull nodes directly (the AbstractNode dispatch) and verify the inside/outside
    # result, including the collinear boundary behaviour.
    coords = [(i, j) for i in 1:5 for j in 1:5]
    for collinear in (false, true)
        @testset "collinear=$collinear" begin
            h = MutableConvexHull{eltype(coords)}(; orientation=CCW, collinear)
            mergepoints!(h, coords)
            # Strictly interior and exterior points are unaffected by the collinear flag.
            @test insidehull((3, 3), h) == true
            @test insidehull((10, 10), h) == false
            # h.hull.head.next is the bottom-left corner (1,1), on the hull boundary.
            # With collinear=false the ≥ test includes the boundary; with collinear=true the > does not.
            @test insidehull(h.hull.head.next, h) == !collinear
        end
    end
end

@testset "Base.in delegates to insidehull" begin
    coords = [(i, j) for i in 1:5 for j in 1:5]
    queries = [(i, j) for i in 0:6 for j in 0:6]
    # ∈ / in must agree with insidehull for every hull kind and collinear setting.
    for H in (MutableConvexHull, MutableLowerConvexHull, MutableUpperConvexHull)
        for collinear in (false, true)
            @testset "$H collinear=$collinear" begin
                h = H{eltype(coords)}(; orientation=CCW, collinear)
                mergepoints!(h, coords)
                @test all((q in h) == insidehull(q, h) for q in queries)
                @test all((q ∈ h) == insidehull(q, h) for q in queries)
            end
        end
    end
    # The full hull of the 1:5 grid is the 1..5 square. Strictly interior points
    # are inside and exterior points are outside regardless of `collinear`; with
    # the default `collinear=false`, boundary points (edges and corners) are inside.
    @testset "interior / exterior" begin
        for collinear in (false, true)
            h = MutableConvexHull{eltype(coords)}(; orientation=CCW, collinear)
            mergepoints!(h, coords)
            @test (3, 3) in h        # strictly interior
            @test !((10, 10) in h)   # exterior
        end
        h = MutableConvexHull{eltype(coords)}(; orientation=CCW, collinear=false)
        mergepoints!(h, coords)
        @test (1, 1) in h            # corner, boundary-inclusive when collinear=false
        @test (1, 3) in h            # edge
    end
end

@testset "removepoint! by value, shared sortedby key" begin
    # With sortedby = first coordinate, (1,0) and (1,2) map to the same key.
    # findpointnode must scan back past one to locate the other.
    sb = x -> x[1]
    pts = [(0, 1), (2, 1), (1, 0), (1, 2)]
    for (H, truthfun) in ((MutableConvexHull,      jarvismarch),
                          (MutableLowerConvexHull,  lower_jarvismarch),
                          (MutableUpperConvexHull,  upper_jarvismarch))
        @testset "$H" begin
            h = H{eltype(pts), typeof(sb)}(; orientation=CCW, collinear=false, sortedby=sb)
            for p in pts
                addpoint!(h, p)
            end
            removepoint!(h, (1, 0))
            @test h == truthfun(filter(!=((1, 0)), pts); sortedby=sb)
        end
    end
end
