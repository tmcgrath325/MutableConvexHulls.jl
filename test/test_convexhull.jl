@testset "convex hulls" begin
    boxcoords = [(i,j) for i in 1:3 for j in 1:3]
    by = x -> (x[1], -x[2])
    n = 10
    
    @testset "convex hull" begin
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
                h1 = MutableLowerConvexHull{eltype(shuffledcoords)}()
                h2 = MutableLowerConvexHull{eltype(shuffledcoords)}(CCW, false, identity, true)
                h3 = MutableLowerConvexHull{eltype(shuffledcoords)}(CW, false, identity, true)
                h4 = MutableLowerConvexHull{eltype(shuffledcoords)}(CCW, true)
                h5 = MutableLowerConvexHull{eltype(shuffledcoords)}(CW, false, by)
                h6 = MutableLowerConvexHull{eltype(shuffledcoords)}(CW, true, by, true)
                hulls = [h1, h2, h3, h4, h5, h6]
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
    
    @testset "lower convex hull" begin
        @testset "initialize" begin
            h = MutableLowerConvexHull{}
        end
        @testset "iterate" begin
            
        end
        @testset "empty" begin
            
        end
        @testset "add point" begin
            
        end
        @testset "merge points" begin

        end
        @testset "remove point" begin

        end
    end

    @testset "upper convex hull" begin
        @testset "initialize" begin
            
        end
        @testset "iterate" begin
            
        end
        @testset "empty" begin
            
        end
        @testset "add point" begin
            
        end
        @testset "merge points" begin

        end
        @testset "remove point" begin

        end
    end
end