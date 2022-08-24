@testset "convex hulls" begin
    boxcoords = [(i,j) for i in 1:3 for j in 1:3]
    by = x -> (x[1], -x[2])
    n = 10
    
    @testset "lower convex hull" begin
        @testset "initialize" begin
            h = MutableLowerConvexHull{eltype(boxcoords)}()
            @test h.orientation === CCW && h.collinear === false && h.sortedby === identity && h.issorted === false
            hull = PairedLinkedList{eltype(boxcoords)}()
            points = PairedLinkedList{eltype(boxcoords)}()
            addpartner!(hull, points)
            h2 = MutableLowerConvexHull{eltype(boxcoords), typeof(by)}(hull, CW, true, by, false)
        end
        @testset "iterate" begin
            
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
                hulls = [MutableLowerConvexHull{eltype(shuffledcoords)}(o, c, f, s) for o in [CCW,CW] for c in [false,true] for f in [identity, by] for s in [false,true]]
                for (i, coord) in enumerate(shuffledcoords)
                    for h in hulls
                        addpoint!(h, coord)
                        @test h == lower_jarvismarch(shuffledcoords[1:i]; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                    end
                end
            end
        end
        @testset "merge points" begin

        end
        @testset "remove point" begin

        end
    end
    
    @testset "upper convex hull" begin
        @testset "initialize" begin
            h = MutableUpperConvexHull{eltype(boxcoords)}()
            @test h.orientation === CCW && h.collinear === false && h.sortedby === identity && h.issorted === false
            hull = PairedLinkedList{eltype(boxcoords)}()
            points = PairedLinkedList{eltype(boxcoords)}()
            addpartner!(hull, points)
            h2 = MutableUpperConvexHull{eltype(boxcoords), typeof(by)}(hull, CW, true, by, false)
        end
        @testset "iterate" begin
            
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
                hulls = [MutableUpperConvexHull{eltype(shuffledcoords)}(o, c, f, s) for o in [CCW,CW] for c in [false,true] for f in [identity, by] for s in [false,true]]
                for (i, coord) in enumerate(shuffledcoords)
                    for h in hulls
                        addpoint!(h, coord)
                        @test h == upper_jarvismarch(shuffledcoords[1:i]; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                    end
                end
            end
        end
        @testset "merge points" begin

        end
        @testset "remove point" begin

        end
    end

    @testset "convex hull" begin
        @testset "initialize" begin
            h = MutableConvexHull{eltype(boxcoords)}()
            @test h.orientation === CCW && h.collinear === false && h.sortedby === identity && h.issorted === false
            hull = PairedLinkedList{eltype(boxcoords)}()
            points = PairedLinkedList{eltype(boxcoords)}()
            addpartner!(hull, points)
            h2 = MutableConvexHull{eltype(boxcoords), typeof(by)}(hull, CW, true, by, false)
        end
        @testset "iterate" begin
            
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
                hulls = [MutableConvexHull{eltype(shuffledcoords)}(o, c, f, s) for o in [CCW,CW] for c in [false,true] for f in [identity, by] for s in [false,true]]
                for (i, coord) in enumerate(shuffledcoords)
                    for h in hulls
                        addpoint!(h, coord)
                        @test h == jarvismarch(shuffledcoords[1:i]; orientation=h.orientation, collinear=h.collinear, sortedby=h.sortedby)
                    end
                end
            end
        end
        @testset "merge points" begin

        end
        @testset "remove point" begin

        end
    end
end