@testset "convex hulls, unique points" begin
    boxcoords = [(i,j) for i in 1:10 for j in 1:10]
    by = x -> (x[1], -x[2])
    n = 10
    
    @testset "lower convex hull" begin
        @testset "initialize" begin
            h = MutableLowerConvexHull{eltype(boxcoords)}()
            @test h.orientation === CCW && h.collinear === false && h.sortedby === identity
            hull = HullList{eltype(boxcoords)}(;sortedby=by)
            points = PointList{eltype(boxcoords)}(;sortedby=by)
            addtarget!(hull, points)
            h2 = MutableLowerConvexHull{eltype(boxcoords), typeof(by)}(hull, points, CW, true, by)
        end
        @testset "empty" begin
            h = lower_monotonechain(boxcoords; orientation=CW, collinear=true, sortedby=by)
            @test length(h) > 0
            h2 = empty(h)
            @test length(h2) == 0 
            @test length(h) > 0
            @test h.orientation == h2.orientation && h.collinear && h2.collinear && h.sortedby == h2.sortedby
            empty!(h)
            @test length(h) == 0
        end
        @testset "add point" begin
            for j=1:n
                shuffledcoords = shuffle(boxcoords)
                hulls = [MutableLowerConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
                for (i, coord) in enumerate(shuffledcoords)
                    for h in hulls
                        addpoint!(h, coord)
                        @test h == lower_jarvismarch(shuffledcoords[1:i]; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                    end
                end
            end
        end
        @testset "merge points" begin
            for j=1:n
                shuffledcoords = shuffle(boxcoords)
                len = Int(sqrt(length(boxcoords)))
                splitcoords = [shuffledcoords[len*(i-1)+1:len*i] for i=1:len]
                hulls = [MutableLowerConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
                mergedcoords = eltype(boxcoords)[]
                for scoords in splitcoords
                    append!(mergedcoords, scoords)
                    for h in hulls
                        mergepoints!(h, scoords)
                        @test h == lower_jarvismarch(mergedcoords; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                    end
                end
            end
        end
        @testset "remove point" begin
            for j=1:n
                shuffledcoords = shuffle(boxcoords)
                hulls = [MutableLowerConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
                for h in hulls
                    for coord in shuffledcoords
                        addpoint!(h, coord)
                    end
                    @test h == lower_jarvismarch(shuffledcoords; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                end
                for i=1:length(shuffledcoords)
                    removeidx = rand(1:length(shuffledcoords))
                    removeddata = shuffledcoords[removeidx]
                    deleteat!(shuffledcoords, removeidx)
                    for h in hulls
                        removepoint!(h, getfirst(x -> x.data == removeddata, ListNodeIterator(h.hull.target)))
                        @test h == lower_jarvismarch(shuffledcoords; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                    end
                end
            end
        end
    end
    
    @testset "upper convex hull" begin
        @testset "initialize" begin
            h = MutableUpperConvexHull{eltype(boxcoords)}()
            @test h.orientation === CCW && h.collinear === false && h.sortedby === identity
            hull = HullList{eltype(boxcoords)}(;sortedby=by)
            points = PointList{eltype(boxcoords)}(;sortedby=by)
            addtarget!(hull, points)
            h2 = MutableUpperConvexHull{eltype(boxcoords), typeof(by)}(hull, points, CW, true, by)
        end
        @testset "empty" begin
            h = upper_monotonechain(boxcoords; orientation=CW, collinear=true, sortedby=by)
            @test length(h) > 0
            h2 = empty(h)
            @test length(h2) == 0 
            @test length(h) > 0
            @test h.orientation == h2.orientation && h.collinear && h2.collinear && h.sortedby == h2.sortedby
            empty!(h)
            @test length(h) == 0
        end
        @testset "add point" begin
            for j=1:n
                shuffledcoords = shuffle(boxcoords)
                hulls = [MutableUpperConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
                for (i, coord) in enumerate(shuffledcoords)
                    for h in hulls
                        addpoint!(h, coord)
                        @test h == upper_jarvismarch(shuffledcoords[1:i]; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                    end
                end
            end
        end
        @testset "merge points" begin
            for j=1:n
                shuffledcoords = shuffle(boxcoords)
                len = Int(sqrt(length(boxcoords)))
                splitcoords = [shuffledcoords[len*(i-1)+1:len*i] for i=1:len]
                hulls = [MutableUpperConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
                mergedcoords = eltype(boxcoords)[]
                for scoords in splitcoords
                    append!(mergedcoords, scoords)
                    for h in hulls
                        mergepoints!(h, scoords)
                        @test h == upper_jarvismarch(mergedcoords; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                    end
                end
            end
        end
        @testset "remove point" begin
            for j=1:n
                shuffledcoords = shuffle(boxcoords)
                hulls = [MutableUpperConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
                for h in hulls
                    for coord in shuffledcoords
                        addpoint!(h, coord)
                    end
                    @test h == upper_jarvismarch(shuffledcoords; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                end
                for i=1:length(shuffledcoords)
                    removeidx = rand(1:length(shuffledcoords))
                    removeddata = shuffledcoords[removeidx]
                    deleteat!(shuffledcoords, removeidx)
                    for h in hulls
                        removepoint!(h, getfirst(x -> x.data == removeddata, ListNodeIterator(h.hull.target)))
                        @test h == upper_jarvismarch(shuffledcoords; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                    end
                end
            end
        end
    end

    @testset "convex hull" begin
        @testset "initialize" begin
            h = MutableConvexHull{eltype(boxcoords)}()
            @test h.orientation === CCW && h.collinear === false && h.sortedby === identity
            hull = HullList{eltype(boxcoords)}(;sortedby=by)
            points = PointList{eltype(boxcoords)}(;sortedby=by)
            addtarget!(hull, points)
            h2 = MutableConvexHull{eltype(boxcoords), typeof(by)}(hull, points, CW, true, by)
        end
        @testset "empty" begin
            h = monotonechain(boxcoords; orientation=CW, collinear=true, sortedby=by)
            @test length(h) > 0
            h2 = empty(h)
            @test length(h2) == 0 
            @test length(h) > 0
            @test h.orientation == h2.orientation && h.collinear && h2.collinear && h.sortedby == h2.sortedby
            empty!(h)
            @test length(h) == 0
        end
        @testset "add point" begin
            for j=1:n
                shuffledcoords = shuffle(boxcoords)
                hulls = [MutableConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
                for (i, coord) in enumerate(shuffledcoords)
                    for h in hulls
                        addpoint!(h, coord)
                        @test h == jarvismarch(shuffledcoords[1:i]; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                    end
                end
            end
        end
        @testset "merge points" begin
            for j=1:n
                shuffledcoords = shuffle(boxcoords)
                len = Int(sqrt(length(boxcoords)))
                splitcoords = [shuffledcoords[len*(i-1)+1:len*i] for i=1:len]
                hulls = [MutableConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
                mergedcoords = eltype(boxcoords)[]
                for scoords in splitcoords
                    append!(mergedcoords, scoords)
                    for h in hulls
                        mergepoints!(h, scoords)
                        @test h == jarvismarch(mergedcoords; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                    end
                end
            end
        end
        @testset "remove point" begin
            for j=1:n
                shuffledcoords = shuffle(boxcoords)
                hulls = [MutableConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
                for h in hulls
                    for coord in shuffledcoords
                        addpoint!(h, coord)
                    end
                    @test h == jarvismarch(shuffledcoords; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                end
                for i=1:length(shuffledcoords)
                    removeidx = rand(1:length(shuffledcoords))
                    removeddata = shuffledcoords[removeidx]
                    deleteat!(shuffledcoords, removeidx)
                    for h in hulls
                        removepoint!(h, getfirst(x -> x.data == removeddata, ListNodeIterator(h.hull.target)))
                    end
                end
            end
        end
    end
end

@testset "convex hulls, duplicate points" begin
    boxcoords = [(i,j) for i in 1:10 for j in 1:10]
    boxcoords = [boxcoords..., boxcoords..., boxcoords...]
    by = x -> (x[1], -x[2])
    n = 10
    
    @testset "lower convex hull" begin
        @testset "initialize" begin
            h = MutableLowerConvexHull{eltype(boxcoords)}()
            @test h.orientation === CCW && h.collinear === false && h.sortedby === identity
            hull = HullList{eltype(boxcoords)}(;sortedby=by)
            points = PointList{eltype(boxcoords)}(;sortedby=by)
            addtarget!(hull, points)
            h2 = MutableLowerConvexHull{eltype(boxcoords), typeof(by)}(hull, points, CW, true, by)
        end
        @testset "empty" begin
            h = lower_monotonechain(boxcoords; orientation=CW, collinear=true, sortedby=by)
            @test length(h) > 0
            h2 = empty(h)
            @test length(h2) == 0 
            @test length(h) > 0
            @test h.orientation == h2.orientation && h.collinear && h2.collinear && h.sortedby == h2.sortedby
            empty!(h)
            @test length(h) == 0
        end
        @testset "add point" begin
            for j=1:n
                shuffledcoords = shuffle(boxcoords)
                hulls = [MutableLowerConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
                for (i, coord) in enumerate(shuffledcoords)
                    for h in hulls
                        addpoint!(h, coord)
                        @test h == lower_jarvismarch(shuffledcoords[1:i]; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                    end
                end
            end
        end
        @testset "merge points" begin
            for j=1:n
                shuffledcoords = shuffle(boxcoords)
                len = Int(sqrt(length(boxcoords) / 3))
                splitcoords = [shuffledcoords[len*(i-1)+1:len*i] for i=1:2*len]
                hulls = [MutableLowerConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
                mergedcoords = eltype(boxcoords)[]
                for scoords in splitcoords
                    append!(mergedcoords, scoords)
                    for h in hulls
                        mergepoints!(h, scoords)
                        @test h == lower_jarvismarch(mergedcoords; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                    end
                end
            end
        end
        @testset "remove point" begin
            for j=1:n
                shuffledcoords = shuffle(boxcoords)
                hulls = [MutableLowerConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
                for h in hulls
                    for coord in shuffledcoords
                        addpoint!(h, coord)
                    end
                    @test h == lower_jarvismarch(shuffledcoords; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                end
                for i=1:length(shuffledcoords)
                    removeidx = rand(1:length(shuffledcoords))
                    removeddata = shuffledcoords[removeidx]
                    deleteat!(shuffledcoords, removeidx)
                    for h in hulls
                        removepoint!(h, getfirst(x -> x.data == removeddata, ListNodeIterator(h.hull.target)))
                        @test h == lower_jarvismarch(shuffledcoords; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                        if h != lower_jarvismarch(shuffledcoords; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                            @show removeddata
                            @show h.subhulls
                            @show h.hull
                            sleep(5)
                        end
                    end
                end
            end
        end
    end
    
    @testset "upper convex hull" begin
        @testset "initialize" begin
            h = MutableUpperConvexHull{eltype(boxcoords)}()
            @test h.orientation === CCW && h.collinear === false && h.sortedby === identity
            hull = HullList{eltype(boxcoords)}(;sortedby=by)
            points = PointList{eltype(boxcoords)}(;sortedby=by)
            addtarget!(hull, points)
            h2 = MutableUpperConvexHull{eltype(boxcoords), typeof(by)}(hull, points, CW, true, by)
        end
        @testset "empty" begin
            h = upper_monotonechain(boxcoords; orientation=CW, collinear=true, sortedby=by)
            @test length(h) > 0
            h2 = empty(h)
            @test length(h2) == 0 
            @test length(h) > 0
            @test h.orientation == h2.orientation && h.collinear && h2.collinear && h.sortedby == h2.sortedby
            empty!(h)
            @test length(h) == 0
        end
        @testset "add point" begin
            for j=1:n
                shuffledcoords = shuffle(boxcoords)
                hulls = [MutableUpperConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
                for (i, coord) in enumerate(shuffledcoords)
                    for h in hulls
                        addpoint!(h, coord)
                        @test h == upper_jarvismarch(shuffledcoords[1:i]; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                    end
                end
            end
        end
        @testset "merge points" begin
            for j=1:n
                shuffledcoords = shuffle(boxcoords)
                len = Int(sqrt(length(boxcoords) / 3))
                splitcoords = [shuffledcoords[len*(i-1)+1:len*i] for i=1:2*len]
                hulls = [MutableUpperConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
                mergedcoords = eltype(boxcoords)[]
                for scoords in splitcoords
                    append!(mergedcoords, scoords)
                    for h in hulls
                        mergepoints!(h, scoords)
                        @test h == upper_jarvismarch(mergedcoords; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                    end
                end
            end
        end
        @testset "remove point" begin
            for j=1:n
                shuffledcoords = shuffle(boxcoords)
                hulls = [MutableUpperConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
                for h in hulls
                    for coord in shuffledcoords
                        addpoint!(h, coord)
                    end
                    @test h == upper_jarvismarch(shuffledcoords; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                end
                for i=1:length(shuffledcoords)
                    removeidx = rand(1:length(shuffledcoords))
                    removeddata = shuffledcoords[removeidx]
                    deleteat!(shuffledcoords, removeidx)
                    for h in hulls
                        removepoint!(h, getfirst(x -> x.data == removeddata, ListNodeIterator(h.hull.target)))
                        @test h == upper_jarvismarch(shuffledcoords; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                    end
                end
            end
        end
    end

    @testset "convex hull" begin
        @testset "initialize" begin
            h = MutableConvexHull{eltype(boxcoords)}()
            @test h.orientation === CCW && h.collinear === false && h.sortedby === identity
            hull = HullList{eltype(boxcoords)}(;sortedby=by)
            points = PointList{eltype(boxcoords)}(;sortedby=by)
            addtarget!(hull, points)
            h2 = MutableConvexHull{eltype(boxcoords), typeof(by)}(hull, points, CW, true, by)
        end
        @testset "empty" begin
            h = monotonechain(boxcoords; orientation=CW, collinear=true, sortedby=by)
            @test length(h) > 0
            h2 = empty(h)
            @test length(h2) == 0 
            @test length(h) > 0
            @test h.orientation == h2.orientation && h.collinear && h2.collinear && h.sortedby == h2.sortedby
            empty!(h)
            @test length(h) == 0
        end
        @testset "add point" begin
            for j=1:n
                shuffledcoords = shuffle(boxcoords)
                hulls = [MutableConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
                for (i, coord) in enumerate(shuffledcoords)
                    for h in hulls
                        addpoint!(h, coord)
                        @test h == jarvismarch(shuffledcoords[1:i]; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                    end
                end
            end
        end
        @testset "merge points" begin
            for j=1:n
                shuffledcoords = shuffle(boxcoords)
                len = Int(sqrt(length(boxcoords) / 3))
                splitcoords = [shuffledcoords[len*(i-1)+1:len*i] for i=1:2*len]
                hulls = [MutableConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
                mergedcoords = eltype(boxcoords)[]
                for scoords in splitcoords
                    append!(mergedcoords, scoords)
                    for h in hulls
                        mergepoints!(h, scoords)
                        @test h == jarvismarch(mergedcoords; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                    end
                end
            end
        end
        @testset "remove point" begin
            for j=1:n
                shuffledcoords = shuffle(boxcoords)
                hulls = [MutableConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
                for h in hulls
                    for coord in shuffledcoords
                        addpoint!(h, coord)
                    end
                    @test h == jarvismarch(shuffledcoords; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                end
                for i=1:length(shuffledcoords)
                    removeidx = rand(1:length(shuffledcoords))
                    removeddata = shuffledcoords[removeidx]
                    deleteat!(shuffledcoords, removeidx)
                    for h in hulls
                        prevhull = collect(h)
                        removepoint!(h, getfirst(x -> x.data == removeddata, ListNodeIterator(h.hull.target)))
                        jarvishull = jarvismarch(shuffledcoords; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                        @test h == jarvishull
                        if h != jarvishull
                            @show prevhull
                            @show shuffledcoords
                            @show removeddata
                            @show h.points
                            @show h.hull
                            @show h.orientation, h.collinear, h.sortedby
                            sleep(5)
                        end
                    end
                end
            end
        end
    end
end

@testset "convex hulls, random data with duplicates" begin
    by = x -> (x[1], -x[2])
    n = 10
    
    @testset "lower convex hull" begin
        @testset "initialize" begin
            boxcoords = [(randn(),randn())]
            h = MutableLowerConvexHull{eltype(boxcoords)}()
            @test h.orientation === CCW && h.collinear === false && h.sortedby === identity
            hull = HullList{eltype(boxcoords)}(;sortedby=by)
            points = PointList{eltype(boxcoords)}(;sortedby=by)
            addtarget!(hull, points)
            h2 = MutableLowerConvexHull{eltype(boxcoords), typeof(by)}(hull, points, CW, true, by)
        end
        @testset "empty" begin
            boxcoords = [(randn(),randn()) for i in 1:10 for j in 1:10]
            h = lower_monotonechain(boxcoords; orientation=CW, collinear=true, sortedby=by)
            @test length(h) > 0
            h2 = empty(h)
            @test length(h2) == 0 
            @test length(h) > 0
            @test h.orientation == h2.orientation && h.collinear && h2.collinear && h.sortedby == h2.sortedby
            empty!(h)
            @test length(h) == 0
        end
        @testset "add point" begin
            for j=1:n
                boxcoords = [(randn(),randn()) for i in 1:10 for j in 1:10]
                boxcoords = [boxcoords..., boxcoords..., boxcoords...]
                shuffledcoords = shuffle(boxcoords)
                hulls = [MutableLowerConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
                for (i, coord) in enumerate(shuffledcoords)
                    for h in hulls
                        addpoint!(h, coord)
                        @test h == lower_jarvismarch(shuffledcoords[1:i]; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                    end
                end
            end
        end
        @testset "merge points" begin
            for j=1:n
                boxcoords = [(randn(),randn()) for i in 1:10 for j in 1:10]
                boxcoords = [boxcoords..., boxcoords..., boxcoords...]
                shuffledcoords = shuffle(boxcoords)
                len = Int(sqrt(length(boxcoords) / 3))
                splitcoords = [shuffledcoords[len*(i-1)+1:len*i] for i=1:2*len]
                hulls = [MutableLowerConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
                mergedcoords = eltype(boxcoords)[]
                for scoords in splitcoords
                    append!(mergedcoords, scoords)
                    for h in hulls
                        mergepoints!(h, scoords)
                        @test h == lower_jarvismarch(mergedcoords; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                    end
                end
            end
        end
        @testset "remove point" begin
            for j=1:n
                boxcoords = [(randn(),randn()) for i in 1:10 for j in 1:10]
                boxcoords = [boxcoords..., boxcoords..., boxcoords...]
                shuffledcoords = shuffle(boxcoords)
                hulls = [MutableLowerConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
                for h in hulls
                    for coord in shuffledcoords
                        addpoint!(h, coord)
                    end
                    @test h == lower_jarvismarch(shuffledcoords; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                end
                for i=1:length(shuffledcoords)
                    removeidx = rand(1:length(shuffledcoords))
                    removeddata = shuffledcoords[removeidx]
                    deleteat!(shuffledcoords, removeidx)
                    for h in hulls
                        removepoint!(h, getfirst(x -> x.data == removeddata, ListNodeIterator(h.hull.target)))
                        @test h == lower_jarvismarch(shuffledcoords; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                        if h != lower_jarvismarch(shuffledcoords; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                            @show removeddata
                            @show h.subhulls
                            @show h.hull
                            sleep(5)
                        end
                    end
                end
            end
        end
    end
    
    @testset "upper convex hull" begin
        @testset "initialize" begin
            boxcoords = [(randn(),randn())]
            h = MutableUpperConvexHull{eltype(boxcoords)}()
            @test h.orientation === CCW && h.collinear === false && h.sortedby === identity
            hull = HullList{eltype(boxcoords)}(;sortedby=by)
            points = PointList{eltype(boxcoords)}(;sortedby=by)
            addtarget!(hull, points)
            h2 = MutableUpperConvexHull{eltype(boxcoords), typeof(by)}(hull, points, CW, true, by)
        end
        @testset "empty" begin
            boxcoords = [(randn(),randn()) for i in 1:10 for j in 1:10]
            h = upper_monotonechain(boxcoords; orientation=CW, collinear=true, sortedby=by)
            @test length(h) > 0
            h2 = empty(h)
            @test length(h2) == 0 
            @test length(h) > 0
            @test h.orientation == h2.orientation && h.collinear && h2.collinear && h.sortedby == h2.sortedby
            empty!(h)
            @test length(h) == 0
        end
        @testset "add point" begin
            for j=1:n
                boxcoords = [(randn(),randn()) for i in 1:10 for j in 1:10]
                boxcoords = [boxcoords..., boxcoords..., boxcoords...]
                shuffledcoords = shuffle(boxcoords)
                hulls = [MutableUpperConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
                for (i, coord) in enumerate(shuffledcoords)
                    for h in hulls
                        addpoint!(h, coord)
                        @test h == upper_jarvismarch(shuffledcoords[1:i]; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                    end
                end
            end
        end
        @testset "merge points" begin
            for j=1:n
                boxcoords = [(randn(),randn()) for i in 1:10 for j in 1:10]
                boxcoords = [boxcoords..., boxcoords..., boxcoords...]
                shuffledcoords = shuffle(boxcoords)
                len = Int(sqrt(length(boxcoords) / 3))
                splitcoords = [shuffledcoords[len*(i-1)+1:len*i] for i=1:2*len]
                hulls = [MutableUpperConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
                mergedcoords = eltype(boxcoords)[]
                for scoords in splitcoords
                    append!(mergedcoords, scoords)
                    for h in hulls
                        mergepoints!(h, scoords)
                        @test h == upper_jarvismarch(mergedcoords; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                    end
                end
            end
        end
        @testset "remove point" begin
            for j=1:n
                boxcoords = [(randn(),randn()) for i in 1:10 for j in 1:10]
                boxcoords = [boxcoords..., boxcoords..., boxcoords...]
                shuffledcoords = shuffle(boxcoords)
                hulls = [MutableUpperConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
                for h in hulls
                    for coord in shuffledcoords
                        addpoint!(h, coord)
                    end
                    @test h == upper_jarvismarch(shuffledcoords; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                end
                for i=1:length(shuffledcoords)
                    removeidx = rand(1:length(shuffledcoords))
                    removeddata = shuffledcoords[removeidx]
                    deleteat!(shuffledcoords, removeidx)
                    for h in hulls
                        removepoint!(h, getfirst(x -> x.data == removeddata, ListNodeIterator(h.hull.target)))
                        @test h == upper_jarvismarch(shuffledcoords; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                    end
                end
            end
        end
    end

    @testset "convex hull" begin
        @testset "initialize" begin
            boxcoords = [(randn(),randn())]
            h = MutableConvexHull{eltype(boxcoords)}()
            @test h.orientation === CCW && h.collinear === false && h.sortedby === identity
            hull = HullList{eltype(boxcoords)}(;sortedby=by)
            points = PointList{eltype(boxcoords)}(;sortedby=by)
            addtarget!(hull, points)
            h2 = MutableConvexHull{eltype(boxcoords), typeof(by)}(hull, points, CW, true, by)
        end
        @testset "empty" begin
            boxcoords = [(randn(),randn()) for i in 1:10 for j in 1:10]
            h = monotonechain(boxcoords; orientation=CW, collinear=true, sortedby=by)
            @test length(h) > 0
            h2 = empty(h)
            @test length(h2) == 0 
            @test length(h) > 0
            @test h.orientation == h2.orientation && h.collinear && h2.collinear && h.sortedby == h2.sortedby
            empty!(h)
            @test length(h) == 0
        end
        @testset "add point" begin
            for j=1:n
                boxcoords = [(randn(),randn()) for i in 1:10 for j in 1:10]
                boxcoords = [boxcoords..., boxcoords..., boxcoords...]
                shuffledcoords = shuffle(boxcoords)
                hulls = [MutableConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
                for (i, coord) in enumerate(shuffledcoords)
                    for h in hulls
                        addpoint!(h, coord)
                        @test h == jarvismarch(shuffledcoords[1:i]; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                    end
                end
            end
        end
        @testset "merge points" begin
            for j=1:n
                boxcoords = [(randn(),randn()) for i in 1:10 for j in 1:10]
                boxcoords = [boxcoords..., boxcoords..., boxcoords...]
                shuffledcoords = shuffle(boxcoords)
                len = Int(sqrt(length(boxcoords) / 3))
                splitcoords = [shuffledcoords[len*(i-1)+1:len*i] for i=1:2*len]
                hulls = [MutableConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
                mergedcoords = eltype(boxcoords)[]
                for scoords in splitcoords
                    append!(mergedcoords, scoords)
                    for h in hulls
                        mergepoints!(h, scoords)
                        @test h == jarvismarch(mergedcoords; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                    end
                end
            end
        end
        @testset "remove point" begin
            for j=1:n
                boxcoords = [(randn(),randn()) for i in 1:10 for j in 1:10]
                boxcoords = [boxcoords..., boxcoords..., boxcoords...]
                shuffledcoords = shuffle(boxcoords)
                hulls = [MutableConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
                for h in hulls
                    for coord in shuffledcoords
                        addpoint!(h, coord)
                    end
                    @test h == jarvismarch(shuffledcoords; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                end
                for i=1:length(shuffledcoords)
                    removeidx = rand(1:length(shuffledcoords))
                    removeddata = shuffledcoords[removeidx]
                    deleteat!(shuffledcoords, removeidx)
                    for h in hulls
                        prevhull = collect(h)
                        removepoint!(h, getfirst(x -> x.data == removeddata, ListNodeIterator(h.hull.target)))
                        jarvishull = jarvismarch(shuffledcoords; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                        @test h == jarvishull
                        if h != jarvishull
                            @show prevhull
                            @show shuffledcoords
                            @show removeddata
                            @show h.points
                            @show h.hull
                            @show h.orientation, h.collinear, h.sortedby
                            sleep(5)
                        end
                    end
                end
            end
        end
    end
end