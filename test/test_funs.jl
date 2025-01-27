# convex hull tests

function hulltestset(name, n, by, coords, hullfun, truthfun)
    @testset "$name" begin
        inittest(by, coords, hullfun)
        emptytest(by, coords, truthfun)
        addtests(n, by, coords, hullfun, truthfun)
        mergetests(n, by, coords, hullfun, truthfun)
        removetests(n, by, coords, hullfun, truthfun)
    end
end

function inittest(by, coords, hullfun)
    @testset "initialize" begin
        h = hullfun{eltype(coords)}()
        @test h.orientation === CCW && h.collinear === false && h.sortedby === identity
        hull = HullList{eltype(coords)}(;sortedby=by)
        points = PointList{eltype(coords)}(;sortedby=by)
        addtarget!(hull, points)
        h2 = hullfun{eltype(coords), typeof(by)}(hull, points, CW, true, by)
    end
end

function emptytest(by, coords, truthfun)
    @testset "empty" begin
        h = truthfun(coords; orientation=CW, collinear=true, sortedby=by)
        @test length(h) > 0
        h2 = empty(h)
        @test length(h2) == 0 
        @test length(h) > 0
        @test h.orientation == h2.orientation && h.collinear && h2.collinear && h.sortedby == h2.sortedby
        empty!(h)
        @test length(h) == 0
    end
end

function addtests(n, by, coords, hullfun, truthfun)
    @testset "add points" begin
        for j=1:n
            shuffledcoords = shuffle(coords)
            hulls = [hullfun{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
            for (i, coord) in enumerate(shuffledcoords)
                for h in hulls
                    addpoint!(h, coord)
                    @test h == truthfun(shuffledcoords[1:i]; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                end
            end
        end
    end
end

function mergetests(n, by, coords, hullfun, truthfun)
    @testset "merge points" begin
        for j=1:n
            shuffledcoords = shuffle(coords)
            mergesize = Int(ceil((sqrt(length(coords))/10)))
            nummerges = Int(length(coords)/mergesize) # have to choose a size of the coords that is divisible by 10
            splitcoords = [shuffledcoords[mergesize*(i-1)+1:mergesize*i] for i=1:nummerges]
            hulls = [hullfun{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
            mergedcoords = eltype(coords)[]
            for scoords in splitcoords
                append!(mergedcoords, scoords)
                for h in hulls
                    mergepoints!(h, scoords)
                    @test h == truthfun(mergedcoords; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                end
            end
        end
    end
end

function removetests(n, by, coords, hullfun, truthfun)
    @testset "remove points" begin
        for j=1:n
            shuffledcoords = shuffle(coords)
            hulls = [hullfun{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
            for h in hulls
                mergepoints!(h, shuffledcoords)
            end
            for i=1:length(shuffledcoords)
                removeidx = rand(1:length(shuffledcoords))
                removeddata = shuffledcoords[removeidx]
                deleteat!(shuffledcoords, removeidx)
                for h in hulls
                    removepoint!(h, getfirst(x -> x.data == removeddata, ListNodeIterator(h.hull.target)))
                    @test h == truthfun(shuffledcoords; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                end
            end
        end
    end
end


# chan hull tests

function chanhulltestset(name, n, by, coords, hullfun, truthfun)
    @testset "$name" begin
        chaninittest(by, coords, hullfun)
        emptytest(by, coords, truthfun)
        addtests(n, by, coords, hullfun, truthfun)
        mergetests(n, by, coords, hullfun, truthfun)
        fallbackmergetest(n, by, coords, hullfun, truthfun)
        chanremovetests(n, by, coords, hullfun, truthfun)
    end
end

function chaninittest(by, coords, hullfun)
    @testset "initialize" begin
        h = hullfun{eltype(coords)}()
        @test h.orientation === CCW && h.collinear === false && h.sortedby === identity
    end
end

function fallbackmergetest(n, by, coords, hullfun, truthfun)
    @testset "fallback merge points" begin
        for j=1:n
            shuffledcoords = shuffle(coords)
            mergesize = Int(ceil((sqrt(length(coords))/10)))
            nummerges = Int(length(coords)/mergesize) # have to choose a size of the coords that is divisible by 10
            splitcoords = [shuffledcoords[mergesize*(i-1)+1:mergesize*i] for i=1:nummerges]
            hulls = [hullfun{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
            mergedcoords = eltype(coords)[]
            for scoords in splitcoords
                append!(mergedcoords, scoords)
                for h in hulls
                    npoints = sum(length, h.subhulls)
                    while npoints > 3 && length(h.subhulls)^2 < npoints
                        push!(h.subhulls, eltype(h.subhulls)(h.orientation, h.collinear, h.sortedby))
                        if (h.subhulls[1].points.cache !== nothing)
                            h.subhulls[end].points.cache = PairedLinkedLists.SkipListCache{T}()
                        end
                    end
                    smallhull = argmin(x->length(x.points),h.subhulls)
                    mergepoints!(smallhull, scoords)
                    MutableConvexHulls.fallback_merge_hull_lists!(h)
                    @test h == truthfun(mergedcoords; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                end
            end
        end
    end
end

function chanremovetests(n, by, coords, hullfun, truthfun)
    @testset "remove point" begin
        for j=1:n
            shuffledcoords = shuffle(coords)
            hulls = [hullfun{eltype(shuffledcoords)}(o, c, f) for o in [CCW,CW] for c in [false,true] for f in [identity, by]]
            for h in hulls
                mergesize = Int(ceil((sqrt(length(shuffledcoords))/10)))
                nummerges = Int(length(shuffledcoords)/mergesize) # have to choose a size of the coords that is divisible by 10
                for k=1:nummerges
                    mergepoints!(h, shuffledcoords[mergesize*(k-1)+1:mergesize*k])
                end
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
                    @test h == truthfun(shuffledcoords; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                end
            end
        end
    end
end