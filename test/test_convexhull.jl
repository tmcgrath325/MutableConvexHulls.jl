@testset "convex hulls, unique points" begin
    boxcoords = [(i, j) for i in 1:10 for j in 1:10]
    by = x -> (x[1], -x[2])
    n = 10

    hulltestset("lower convex hull", n, by, boxcoords, MutableLowerConvexHull, lower_jarvismarch)
    hulltestset("upper convex hull", n, by, boxcoords, MutableUpperConvexHull, upper_jarvismarch)
    hulltestset("convex hull", n, by, boxcoords, MutableConvexHull, jarvismarch)
end

@testset "convex hulls, duplicate points" begin
    boxcoords = [(i, j) for i in 1:10 for j in 1:10]
    boxcoords = [boxcoords..., boxcoords..., boxcoords...]
    by = x -> (x[1], -x[2])
    n = 10

    hulltestset("lower convex hull", n, by, boxcoords, MutableLowerConvexHull, lower_jarvismarch)
    hulltestset("upper convex hull", n, by, boxcoords, MutableUpperConvexHull, upper_jarvismarch)
    hulltestset("convex hull", n, by, boxcoords, MutableConvexHull, jarvismarch)
end

@testset "convex hulls, random data with duplicates" begin
    coords = [(randn(), randn()) for i in 1:10 for j in 1:10]
    by = x -> (x[1], -x[2])
    n = 10

    hulltestset("lower convex hull", n, by, coords, MutableLowerConvexHull, lower_jarvismarch)
    hulltestset("upper convex hull", n, by, coords, MutableUpperConvexHull, upper_jarvismarch)
    hulltestset("convex hull", n, by, coords, MutableConvexHull, jarvismarch)
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
        a = MutableConvexHull{Tuple{Float64, Float64}}()
        b = MutableConvexHull{Tuple{Float64, Float64}}()
        for p in [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
            addpoint!(a, p)
        end
        for p in [(2.0, 2.0), (3.0, 2.0), (3.0, 3.0), (2.0, 3.0)]
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

    @testset "mergehulls! rejects Chan hulls" begin
        # Chan hulls lack a .points field; they must not silently crash with a
        # field-access error — the caller should get a MethodError instead.
        T = Tuple{Float64, Float64}
        for H in (ChanConvexHull, ChanLowerConvexHull, ChanUpperConvexHull)
            h = H{T}()
            @test_throws MethodError mergehulls!(h, h)
            @test_throws MethodError mergehulls(h, h)
        end
    end

    @testset "merge_hull_lists! emits no fallback warning" begin
        # The optimized merge path must never silently fall back on valid input.
        pts_a = [(Float64(i), Float64(j)) for i in 1:5 for j in 1:5]
        pts_b = [(Float64(i), Float64(j)) for i in 6:10 for j in 1:5]
        for H in (MutableConvexHull, MutableLowerConvexHull, MutableUpperConvexHull)
            a = H{eltype(pts_a)}(); mergepoints!(a, pts_a)
            b = H{eltype(pts_b)}(); mergepoints!(b, pts_b)
            @test_logs min_level = Logging.Warn mergehulls!(copy(a), copy(b))
        end
        for H in (ChanConvexHull, ChanLowerConvexHull, ChanUpperConvexHull)
            h = H{Tuple{Float64, Float64}}()
            mergepoints!(h, [pts_a..., pts_b...])
            @test_logs min_level = Logging.Warn mergepoints!(h, [(11.0, Float64(j)) for j in 1:5])
        end
    end
end

@testset "constructor keyword API" begin
    # All configuration arguments are keyword-only; positional calls now throw.
    T = Tuple{Float64, Float64}
    by = x -> x[1]
    for H in (MutableConvexHull, MutableLowerConvexHull, MutableUpperConvexHull)
        # Defaults
        h = H{T}()
        @test h.orientation === CCW && h.collinear === false && h.sortedby === identity
        # All keywords explicit
        h2 = H{T}(; orientation = CW, collinear = true, sortedby = by)
        @test h2.orientation === CW && h2.collinear === true && h2.sortedby === by
        # Positional call now throws
        @test_throws MethodError H{T}(CCW, false, identity)
    end
    for H in (ChanConvexHull, ChanLowerConvexHull, ChanUpperConvexHull)
        h = H{T}()
        @test h.orientation === CCW && h.collinear === false && h.sortedby === identity
        h2 = H{T}(; orientation = CW, collinear = true, sortedby = by)
        @test h2.orientation === CW && h2.collinear === true && h2.sortedby === by
        @test_throws MethodError H{T}(CCW, false, identity)
    end
end

@testset "removepoint! by value" begin
    boxcoords = [(i, j) for i in 1:10 for j in 1:10]
    coords = [boxcoords..., boxcoords...]   # include duplicate points
    for (H, truthfun) in (
            (MutableLowerConvexHull, lower_jarvismarch),
            (MutableUpperConvexHull, upper_jarvismarch),
            (MutableConvexHull, jarvismarch),
        )
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
    for (H, truthfun) in (
            (MutableLowerConvexHull, lower_jarvismarch),
            (MutableUpperConvexHull, upper_jarvismarch),
            (MutableConvexHull, jarvismarch),
        )
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
            h = MutableConvexHull{eltype(coords)}(; orientation = CCW, collinear)
            mergepoints!(h, coords)
            # Strictly interior and exterior points are unaffected by the collinear flag.
            @test insidehull((3, 3), h) == true
            @test insidehull((10, 10), h) == false
            # h.hull.head.next is the bottom-left corner (1,1), an exact hull vertex.
            # Exact vertex matches are always inside regardless of collinear.
            @test insidehull(h.hull.head.next, h)
        end
    end
end

@testset "insidehull interior / exterior / boundary" begin
    coords = [(i, j) for i in 1:5 for j in 1:5]
    # With collinear=false the hull is the 4-corner square; with collinear=true it
    # includes all collinear boundary grid points as explicit vertices.
    for collinear in (false, true)
        h = MutableConvexHull{eltype(coords)}(; orientation = CCW, collinear)
        mergepoints!(h, coords)
        @test insidehull((3, 3), h)                 # strictly interior
        @test !insidehull((10, 10), h)              # exterior
        @test insidehull((1, 1), h)                 # corner vertex: always inside regardless of collinear
        # (1,3) is in the input data: a hull vertex when collinear=true, an edge point when
        # collinear=false.  In both cases insidehull returns true.
        @test insidehull((1, 3), h)
    end

    # Non-vertex boundary points (not in the input data): inside when collinear=false,
    # outside when collinear=true.
    fcoords = [(Float64(i), Float64(j)) for i in 1:5 for j in 1:5]
    for collinear in (false, true)
        hf = MutableConvexHull{eltype(fcoords)}(; orientation = CCW, collinear)
        mergepoints!(hf, fcoords)
        @test insidehull((1.0, 2.5), hf) != collinear    # left edge, between (1,2) and (1,3)
    end
end

@testset "insidehull collinear: interior hull vertices always inside" begin
    # Hull vertices that are not on vertical extreme edges previously returned false
    # with collinear=true.  All hull vertices must return true regardless of collinear.
    pts = [(0.0, 1.0), (2.0, 0.0), (4.0, 1.0), (3.0, 3.0), (1.0, 3.0)]
    for H in (
            MutableConvexHull, MutableLowerConvexHull, MutableUpperConvexHull,
            ChanConvexHull, ChanLowerConvexHull, ChanUpperConvexHull,
        )
        for collinear in (false, true)
            h = H{eltype(pts)}(; collinear)
            for p in pts
                addpoint!(h, p)
            end
            for v in collect(h)
                @test insidehull(v, h)
            end
        end
    end
end

@testset "removepoint! by value, shared sortedby key" begin
    # With sortedby = first coordinate, (1,0) and (1,2) map to the same key.
    # findpointnode must scan back past one to locate the other.
    sb = x -> x[1]
    pts = [(0, 1), (2, 1), (1, 0), (1, 2)]
    for (H, truthfun) in (
            (MutableConvexHull, jarvismarch),
            (MutableLowerConvexHull, lower_jarvismarch),
            (MutableUpperConvexHull, upper_jarvismarch),
        )
        @testset "$H" begin
            h = H{eltype(pts), typeof(sb)}(; orientation = CCW, collinear = false, sortedby = sb)
            for p in pts
                addpoint!(h, p)
            end
            removepoint!(h, (1, 0))
            @test h == truthfun(filter(!=((1, 0)), pts); sortedby = sb)
        end
    end
end
