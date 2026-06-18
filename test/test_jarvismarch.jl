# using MutableConvexHulls: jarvismarch!, lower_jarvismarch!, upper_jarvismarch!

@testset "Jarvis March Algorithm" begin
    irange = 1:10
    jrange = 1:10
    boxcoords = [(i, j) for i in irange for j in jrange]
    dupcoords = [boxcoords..., boxcoords...]
    by = x -> (x[1], -x[2])

    @testset "Lower Jarvis March" begin
        # standard-sorted coords
        lower = lower_jarvismarch(boxcoords) # ; orientation = CCW, collinear = false
        lowerCW = lower_jarvismarch(boxcoords; orientation = CW)
        @test collect(lower.hull) == [first(boxcoords), (last(boxcoords)[1], first(boxcoords)[2]), last(boxcoords)]
        @test collect(lowerCW.hull) == reverse(collect(lower.hull))
        lowercollinear = lower_jarvismarch(boxcoords; collinear = true)
        lowercollinearCW = lower_jarvismarch(boxcoords; orientation = CW, collinear = true)
        @test collect(lowercollinear.hull) == [[(i, first(jrange)) for i in irange]..., [(last(irange), j) for j in jrange[2:end]]...]
        @test collect(lowercollinearCW.hull) == reverse(collect(lowercollinear.hull))
        for l in [lower, lowerCW, lowercollinear, lowercollinearCW]
            @test collect(l.hull) == collect(MCH.jarvismarch!(l).hull)
        end

        # alternate-sorting coords
        lower2 = lower_jarvismarch(boxcoords; sortedby = by) # ; orientation = CCW, collinear = false
        lowerCW2 = lower_jarvismarch(boxcoords; orientation = CW, sortedby = by)
        @test collect(lower2.hull) == [(first(boxcoords)[1], last(boxcoords)[2]), first(boxcoords), (last(boxcoords)[1], first(boxcoords)[2])]
        @test collect(lowerCW2.hull) == reverse(collect(lower2.hull))
        lowercollinear2 = lower_jarvismarch(boxcoords; collinear = true, sortedby = by)
        lowercollinearCW2 = lower_jarvismarch(boxcoords; orientation = CW, collinear = true, sortedby = by)
        @test collect(lowercollinear2.hull) == [[(first(jrange), i) for i in reverse(irange)]..., [(j, first(irange)) for j in jrange[2:end]]...]
        @test collect(lowercollinearCW2.hull) == reverse(collect(lowercollinear2.hull))
        for l in [lower2, lowerCW2, lowercollinear2, lowercollinearCW2]
            @test collect(l.hull) == collect(MCH.jarvismarch!(l).hull)
        end
    end

    @testset "Upper Jarvis March" begin
        # standard-sorted coords
        upper = upper_jarvismarch(boxcoords) # ; orientation = CCW, collinear = false
        upperCW = upper_jarvismarch(boxcoords; orientation = CW)
        @test collect(upper.hull) == [last(boxcoords), (first(boxcoords)[1], last(boxcoords)[2]), first(boxcoords)]
        @test collect(upperCW.hull) == reverse(collect(upper.hull))
        uppercollinear = upper_jarvismarch(boxcoords; collinear = true)
        uppercollinearCW = upper_jarvismarch(boxcoords; orientation = CW, collinear = true)
        @test collect(uppercollinear.hull) == [[(i, last(jrange)) for i in reverse(irange)]..., [(first(irange), j) for j in reverse(jrange)[2:end]]...]
        @test collect(uppercollinearCW.hull) == reverse(collect(uppercollinear.hull))
        for u in [upper, upperCW, uppercollinear, uppercollinearCW]
            @test collect(u.hull) == collect(MCH.jarvismarch!(u).hull)
        end

        # alternate-sorting coords
        upper2 = upper_jarvismarch(boxcoords; sortedby = by) # ; orientation = CCW, collinear = false
        upperCW2 = upper_jarvismarch(boxcoords; orientation = CW, sortedby = by)
        @test collect(upper2.hull) == [(last(boxcoords)[1], first(boxcoords)[2]), last(boxcoords), (first(boxcoords)[1], last(boxcoords)[2])]
        @test collect(upperCW2.hull) == reverse(collect(upper2.hull))
        uppercollinear2 = upper_jarvismarch(boxcoords; collinear = true, sortedby = by)
        uppercollinearCW2 = upper_jarvismarch(boxcoords; orientation = CW, collinear = true, sortedby = by)
        @test collect(uppercollinear2.hull) == [[(last(jrange), i) for i in irange]..., [(j, last(irange)) for j in reverse(jrange)[2:end]]...]
        @test collect(uppercollinearCW2.hull) == reverse(collect(uppercollinear2.hull))
        for u in [upper2, upperCW2, uppercollinear2, uppercollinearCW2]
            @test collect(u.hull) == collect(MCH.jarvismarch!(u).hull)
        end
    end

    @testset "Full Jarvis March" begin
        # standard-sorted coords
        hull = jarvismarch(boxcoords) # ; orientation = CCW, collinear = false
        hullCW = jarvismarch(boxcoords; orientation = CW)
        @test collect(hull.hull) == [first(boxcoords), (last(boxcoords)[1], first(boxcoords)[2]), last(boxcoords), (first(boxcoords)[1], last(boxcoords)[2])]
        @test collect(hullCW.hull) == reverse(circshift(collect(hull.hull), Int(length(hull.hull) / 2 - 1)))
        hullcollinear = jarvismarch(boxcoords; collinear = true)
        hullcollinearCW = jarvismarch(boxcoords; orientation = CW, collinear = true)
        @test collect(hullcollinear.hull) == [[(i, first(jrange)) for i in irange]..., [(last(irange), j) for j in jrange[2:end]]..., [(i, last(jrange)) for i in reverse(irange)[2:end]]..., [(first(irange), j) for j in reverse(jrange)[2:(end - 1)]]...]
        @test collect(hullcollinearCW.hull) == reverse(circshift(collect(hullcollinear.hull), Int(length(hullcollinear.hull) / 2) - 1))
        for h in [hull, hullCW, hullcollinear, hullcollinearCW]
            @test collect(h.hull) == collect(MCH.jarvismarch!(h).hull)
        end

        # alternate-sorting coords
        hull2 = jarvismarch(boxcoords; sortedby = by) # ; orientation = CCW, collinear = false
        hullCW2 = jarvismarch(boxcoords; orientation = CW, sortedby = by)
        @test collect(hull2.hull) == [(first(boxcoords)[1], last(boxcoords)[2]), first(boxcoords), (last(boxcoords)[1], first(boxcoords)[2]), last(boxcoords)]
        @test collect(hullCW2.hull) == reverse(circshift(collect(hull2.hull), Int(length(hull2.hull) / 2 - 1)))
        hullcollinear2 = jarvismarch(boxcoords; collinear = true, sortedby = by)
        hullcollinearCW2 = jarvismarch(boxcoords; orientation = CW, collinear = true, sortedby = by)
        @test collect(hullcollinear2.hull) == [[(first(jrange), i) for i in reverse(irange)]..., [(j, first(irange)) for j in jrange[2:end]]..., [(last(jrange), i) for i in irange[2:end]]..., [(j, last(irange)) for j in reverse(jrange)[2:(end - 1)]]...]
        @test collect(hullcollinearCW2.hull) == reverse(circshift(collect(hullcollinear2.hull), Int(length(hullcollinear2.hull) / 2 - 1)))
        for h in [hull2, hullCW2, hullcollinear2, hullcollinearCW2]
            @test collect(h.hull) == collect(MCH.jarvismarch!(h).hull)
        end
    end

    @testset "Matrix input" begin
        m = [p[k] for p in boxcoords, k in 1:2]
        @test lower_jarvismarch(m) == lower_jarvismarch(boxcoords)
        @test upper_jarvismarch(m) == upper_jarvismarch(boxcoords)
        @test jarvismarch(m) == jarvismarch(boxcoords)

        # configuration keywords are forwarded through the matrix form
        @test collect(jarvismarch(m; orientation = CW).hull) == collect(jarvismarch(boxcoords; orientation = CW).hull)
        @test collect(lower_jarvismarch(m; collinear = true).hull) == collect(lower_jarvismarch(boxcoords; collinear = true).hull)
        @test collect(upper_jarvismarch(m; sortedby = by).hull) == collect(upper_jarvismarch(boxcoords; sortedby = by).hull)
    end

    @testset "jarvismarch! return type — 0 and 1 point" begin
        T = Tuple{Int, Int}
        for H in (MutableConvexHull, MutableLowerConvexHull, MutableUpperConvexHull)
            h0 = H{T}()
            @test MCH.jarvismarch!(h0) isa H
            h1 = H{T}()
            addpoint!(h1, (1, 1))
            @test MCH.jarvismarch!(h1) isa H
        end
    end

    @testset "jarvismarch! rejects an empty hull" begin
        T = Tuple{Int, Int}
        for H in (MutableConvexHull, MutableLowerConvexHull, MutableUpperConvexHull)
            h = H{T}()
            @test_throws ArgumentError MCH.jarvismarch!(h.hull, h.hull.target, h.collinear, h.orientation, MCH.DOWN)
            @test_throws "at least one point" MCH.jarvismarch!(h.hull, h.hull.target, h.collinear, h.orientation, MCH.DOWN)
        end
    end
end

@testset "jarvissortedsearch" begin
    # jarvissortedsearch requires a strictly-convex partial (lower or upper) hull list
    # built with collinear=false. Only that domain guarantees the unimodality of the
    # betterturn quality function that the early-exit relies on.
    boxcoords = [(i, j) for i in 1:10 for j in 1:10]
    lhalf     = [(i, j) for i in 1:5  for j in 1:10]
    rhalf     = [(i, j) for i in 6:10 for j in 1:10]

    for orientation in [CCW, CW]
        bt = let orientation = orientation
            (pe, o, a, b) -> MCH.isfurtherturn(!orientation, pe, o, a, b)
        end
        # Use hull nodes from h1 as queries; search within h2's sorted hull list.
        for (coords1, coords2) in [(lhalf, rhalf), (rhalf, lhalf), (lhalf, boxcoords), (boxcoords, rhalf)]
            for H in (MutableLowerConvexHull, MutableUpperConvexHull)
                h1 = H{eltype(coords1)}(; orientation, collinear = false)
                h2 = H{eltype(coords2)}(; orientation, collinear = false)
                mergepoints!(h1, coords1)
                mergepoints!(h2, coords2)
                (isempty(h1.hull) || isempty(h2.hull)) && continue
                for query in ListNodeIterator(h1.hull)
                    for pe in [MCH.UP, MCH.DOWN, (1.0, 0.0), (-1.0, 0.0)]
                        ref    = MCH.jarvissearch(query, pe, ListNodeIterator(h2.hull), bt)
                        result = MCH.jarvissortedsearch(query, pe, h2.hull, bt)
                        @test MCH.coordsareequal(result.data, ref.data)
                    end
                end
            end
        end
    end

    @testset "edge cases" begin
        T = Tuple{Int, Int}
        bt_base = (pe, o, a, b) -> MCH.isfurtherturn(!CCW, pe, o, a, b)

        # Build a one-point hull to get a query node; its hull list provides the query.
        q_hull = MutableLowerConvexHull{T}()
        addpoint!(q_hull, (0, 0))
        query = head(q_hull.hull)

        # Empty list throws.
        @test_throws ArgumentError MCH.jarvissortedsearch(query, MCH.DOWN, MutableLowerConvexHull{T}().hull, bt_base)
        @test_throws "must not be empty" MCH.jarvissortedsearch(query, MCH.DOWN, MutableLowerConvexHull{T}().hull, bt_base)

        # Single-point list returns that point regardless of query or prevedge.
        one_hull = MutableLowerConvexHull{T}()
        addpoint!(one_hull, (5, 3))
        for pe in [MCH.UP, MCH.DOWN, (1.0, 0.0)]
            @test MCH.coordsareequal(MCH.jarvissortedsearch(query, pe, one_hull.hull, bt_base).data, (5, 3))
        end

        # Two-point list: result must match jarvissearch.
        two_hull = MutableLowerConvexHull{T}()
        mergepoints!(two_hull, [(2, 0), (8, 0)])
        for pe in [MCH.UP, MCH.DOWN, (1.0, 0.0), (-1.0, 0.0)]
            ref    = MCH.jarvissearch(query, pe, ListNodeIterator(two_hull.hull), bt_base)
            result = MCH.jarvissortedsearch(query, pe, two_hull.hull, bt_base)
            @test MCH.coordsareequal(result.data, ref.data)
        end

        # List where all points share the query's coordinates: returns head.
        dup_hull = MutableLowerConvexHull{T}()
        addpoint!(dup_hull, (0, 0))  # same coords as query
        @test MCH.coordsareequal(MCH.jarvissortedsearch(query, MCH.DOWN, dup_hull.hull, bt_base).data, (0, 0))
    end
end
