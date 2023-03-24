@testset "chan hulls, unique points" begin
    boxcoords = [(i,j) for i in 1:10 for j in 1:10]
    by = x -> (x[1], -x[2])
    n = 10
    
    @testset "lower chan hull" begin
        @testset "initialize" begin
            h = ChanLowerConvexHull{eltype(boxcoords)}()
            @test h.orientation === CCW && h.collinear === false && h.sortedby === identity
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
                hulls = [ChanLowerConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
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
                hulls = [ChanLowerConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
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
                hulls = [ChanLowerConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
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
                        rmnode = first(h.subhulls).points.head
                        for sh in h.subhulls
                            rmnode = getfirst(x -> x.data == removeddata, ListNodeIterator(sh.points))
                            !isnothing(rmnode) && break
                        end
                        removepoint!(h, rmnode)
                        @test h == lower_jarvismarch(shuffledcoords; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                        if h != lower_jarvismarch(shuffledcoords; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                            @show shuffledcoords
                            @show removeddata
                            @show h.orientation, h.collinear, h.sortedby
                            sleep(5)
                        end
                    end
                end
            end
        end
    end
    
    @testset "upper chan hull" begin
        @testset "initialize" begin
            h = ChanUpperConvexHull{eltype(boxcoords)}()
            @test h.orientation === CCW && h.collinear === false && h.sortedby === identity
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
                hulls = [ChanUpperConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
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
                hulls = [ChanUpperConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
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
                hulls = [ChanUpperConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
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
                        rmnode = first(h.subhulls).points.head
                        for sh in h.subhulls
                            rmnode = getfirst(x -> x.data == removeddata, ListNodeIterator(sh.points))
                            !isnothing(rmnode) && break
                        end
                        removepoint!(h, rmnode)
                        @test h == upper_jarvismarch(shuffledcoords; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                    end
                end
            end
        end
    end

    @testset "chan hull" begin
        @testset "initialize" begin
            h = ChanConvexHull{eltype(boxcoords)}()
            @test h.orientation === CCW && h.collinear === false && h.sortedby === identity
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
                hulls = [ChanConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
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
                hulls = [ChanConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
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
            for j=1:10
                shuffledcoords = shuffle(boxcoords)
                hulls = [ChanConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
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
                        rmnode = first(h.subhulls).points.head
                        for sh in h.subhulls
                            rmnode = getfirst(x -> x.data == removeddata, ListNodeIterator(sh.points))
                            !isnothing(rmnode) && break
                        end
                        removepoint!(h, rmnode)
                        @test h == jarvismarch(shuffledcoords; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                    end
                end
            end
        end
    end
end

@testset "chan hulls, duplicate points" begin
    boxcoords = [(i,j) for i in 1:10 for j in 1:10]
    boxcoords = [boxcoords..., boxcoords..., boxcoords...]
    by = x -> (x[1], -x[2])
    n = 10
    
    @testset "lower chan hull" begin
        @testset "initialize" begin
            h = ChanLowerConvexHull{eltype(boxcoords)}()
            @test h.orientation === CCW && h.collinear === false && h.sortedby === identity
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
                hulls = [ChanLowerConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
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
                hulls = [ChanLowerConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
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
                hulls = [ChanLowerConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
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
                        rmnode = first(h.subhulls).points.head
                        for sh in h.subhulls
                            rmnode = getfirst(x -> x.data == removeddata, ListNodeIterator(sh.points))
                            !isnothing(rmnode) && break
                        end
                        removepoint!(h, rmnode)
                        @test h == lower_jarvismarch(shuffledcoords; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                    end
                end
            end
        end
    end
    
    @testset "upper chan hull" begin
        @testset "initialize" begin
            h = ChanUpperConvexHull{eltype(boxcoords)}()
            @test h.orientation === CCW && h.collinear === false && h.sortedby === identity
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
                hulls = [ChanUpperConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
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
                hulls = [ChanUpperConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
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
                hulls = [ChanUpperConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
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
                        rmnode = first(h.subhulls).points.head
                        for sh in h.subhulls
                            rmnode = getfirst(x -> x.data == removeddata, ListNodeIterator(sh.points))
                            !isnothing(rmnode) && break
                        end
                        removepoint!(h, rmnode)
                        @test h == upper_jarvismarch(shuffledcoords; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                    end
                end
            end
        end
    end

    @testset "chan hull" begin
        @testset "initialize" begin
            h = ChanConvexHull{eltype(boxcoords)}()
            @test h.orientation === CCW && h.collinear === false && h.sortedby === identity
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
                hulls = [ChanConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
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
                hulls = [ChanConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
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
            for j=1:10
                shuffledcoords = shuffle(boxcoords)
                hulls = [ChanConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
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
                        rmnode = first(h.subhulls).points.head
                        for sh in h.subhulls
                            rmnode = getfirst(x -> x.data == removeddata, ListNodeIterator(sh.points))
                            !isnothing(rmnode) && break
                        end
                        removepoint!(h, rmnode)
                        @test h == jarvismarch(shuffledcoords; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                    end
                end
            end
        end
    end
end

@testset "chan hulls, random data with duplicates" begin
    by = x -> (x[1], -x[2])
    n = 10
    
    @testset "lower chan hull" begin
        @testset "initialize" begin
            boxcoords = [(randn(),randn())]
            h = ChanLowerConvexHull{eltype(boxcoords)}()
            @test h.orientation === CCW && h.collinear === false && h.sortedby === identity
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
                hulls = [ChanLowerConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
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
                hulls = [ChanLowerConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
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
                hulls = [ChanLowerConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
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
                        rmnode = first(h.subhulls).points.head
                        for sh in h.subhulls
                            rmnode = getfirst(x -> x.data == removeddata, ListNodeIterator(sh.points))
                            !isnothing(rmnode) && break
                        end
                        removepoint!(h, rmnode)
                        @test h == lower_jarvismarch(shuffledcoords; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                    end
                end
            end
        end
    end
    
    @testset "upper chan hull" begin
        @testset "initialize" begin
            boxcoords = [(randn(),randn())]
            h = ChanUpperConvexHull{eltype(boxcoords)}()
            @test h.orientation === CCW && h.collinear === false && h.sortedby === identity
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
                hulls = [ChanUpperConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
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
                hulls = [ChanUpperConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
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
                hulls = [ChanUpperConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
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
                        rmnode = first(h.subhulls).points.head
                        for sh in h.subhulls
                            rmnode = getfirst(x -> x.data == removeddata, ListNodeIterator(sh.points))
                            !isnothing(rmnode) && break
                        end
                        removepoint!(h, rmnode)
                        @test h == upper_jarvismarch(shuffledcoords; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                    end
                end
            end
        end
    end

    @testset "chan hull" begin
        @testset "initialize" begin
            boxcoords = [(randn(),randn())]
            h = ChanConvexHull{eltype(boxcoords)}()
            @test h.orientation === CCW && h.collinear === false && h.sortedby === identity
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
                hulls = [ChanConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
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
                hulls = [ChanConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
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
                hulls = [ChanConvexHull{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
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
                        rmnode = first(h.subhulls).points.head
                        for sh in h.subhulls
                            rmnode = getfirst(x -> x.data == removeddata, ListNodeIterator(sh.points))
                            !isnothing(rmnode) && break
                        end
                        removepoint!(h, rmnode)
                        @test h == jarvismarch(shuffledcoords; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                    end
                end
            end
        end
    end
end