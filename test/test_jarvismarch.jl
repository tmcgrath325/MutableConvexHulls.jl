# using MutableConvexHulls: jarvismarch!, lower_jarvismarch!, upper_jarvismarch!

@testset "Jarvis March Algorithm" begin
    irange = 1:10
    jrange = 1:10
    boxcoords = [(i,j) for i in irange for j in jrange]
    by = x -> (x[1], -x[2])

    @testset "Lower Jarvis March" begin
        # standard-sorted coords
        lower = lower_jarvismarch(boxcoords) # ; orientation = CCW, collinear = false
        @test collect(lower.hull) == [first(boxcoords), (last(boxcoords)[1], first(boxcoords)[2]), last(boxcoords)]
        @test collect(lower_jarvismarch(boxcoords; orientation=CW).hull) == reverse(collect(lower.hull))
        lowercollinear = lower_jarvismarch(boxcoords; collinear = true)
        @test collect(lowercollinear.hull) == [[(i,first(jrange)) for i in irange]..., [(last(irange), j) for j in jrange[2:end]]...]
        @test collect(lower_jarvismarch(boxcoords; orientation=CW, collinear=true).hull) == reverse(collect(lowercollinear.hull))

        # alternate-sorting coords
        lower2 = lower_jarvismarch(boxcoords; by=by) # ; orientation = CCW, collinear = false
        @test collect(lower2.hull) == [(first(boxcoords)[1], last(boxcoords)[2]), first(boxcoords), (last(boxcoords)[1], first(boxcoords)[2])] 
        @test collect(lower_jarvismarch(boxcoords; orientation=CW, by=by).hull) == reverse(collect(lower2.hull))
        lowercollinear2 = lower_jarvismarch(boxcoords; collinear=true, by=by)
        @test collect(lowercollinear2.hull) == [[(first(jrange),i) for i in reverse(irange)]..., [(j,first(irange)) for j in jrange[2:end]]...] 
        @test collect(lower_jarvismarch(boxcoords; orientation=CW, collinear=true, by=by).hull) == reverse(collect(lowercollinear2.hull))
    end

    @testset "Upper Jarvis March" begin
        # standard-sorted coords
        upper = upper_jarvismarch(boxcoords) # ; orientation = CCW, collinear = false
        @test collect(upper.hull) == [last(boxcoords), (first(boxcoords)[1], last(boxcoords)[2]), first(boxcoords)]
        @test collect(upper_jarvismarch(boxcoords; orientation=CW).hull) == reverse(collect(upper.hull))
        uppercollinear = upper_jarvismarch(boxcoords; collinear = true)
        @test collect(uppercollinear.hull) == [[(i,last(jrange)) for i in reverse(irange)]..., [(first(irange), j) for j in reverse(jrange)[2:end]]...]
        @test collect(upper_jarvismarch(boxcoords; orientation=CW, collinear=true).hull) == reverse(collect(uppercollinear.hull))

        # alternate-sorting coords
        upper2 = upper_jarvismarch(boxcoords; by=by) # ; orientation = CCW, collinear = false
        @test collect(upper2.hull) == [(last(boxcoords)[1], first(boxcoords)[2]), last(boxcoords), (first(boxcoords)[1], last(boxcoords)[2])]
        @test collect(upper_jarvismarch(boxcoords; orientation=CW, by=by).hull) == reverse(collect(upper2.hull))
        uppercollinear2 = upper_jarvismarch(boxcoords; collinear=true, by=by)
        @test collect(uppercollinear2.hull) == [[(last(jrange),i) for i in irange]..., [(j,last(irange)) for j in reverse(jrange)[2:end]]...]
        @test collect(upper_jarvismarch(boxcoords; orientation=CW, collinear=true, by=by).hull) == reverse(collect(uppercollinear2.hull))
    end

    @testset "Full Jarvis March" begin
        # standard-sorted coords
        hull = jarvismarch(boxcoords) # ; orientation = CCW, collinear = false
        @test collect(hull.hull) == [first(boxcoords), (last(boxcoords)[1], first(boxcoords)[2]), last(boxcoords), (first(boxcoords)[1], last(boxcoords)[2])]
        @test collect(jarvismarch(boxcoords; orientation=CW).hull) == reverse(circshift(collect(hull.hull),Int(length(hull.hull)/2-1))) 
        hullcollinear = jarvismarch(boxcoords; collinear = true)
        @test collect(hullcollinear.hull) == [[(i,first(jrange)) for i in irange]..., [(last(irange), j) for j in jrange[2:end]]..., [(i,last(jrange)) for i in reverse(irange)[2:end]]..., [(first(irange), j) for j in reverse(jrange)[2:end-1]]...]
        @test collect(jarvismarch(boxcoords; orientation=CW, collinear=true).hull) == reverse(circshift(collect(hullcollinear.hull), Int(length(hullcollinear.hull)/2)-1))

        # alternate-sorting coords
        hull2 = jarvismarch(boxcoords; by=by) # ; orientation = CCW, collinear = false
        @test collect(hull2.hull) == [(first(boxcoords)[1], last(boxcoords)[2]), first(boxcoords), (last(boxcoords)[1], first(boxcoords)[2]), last(boxcoords)]
        @test collect(jarvismarch(boxcoords; orientation=CW, by=by).hull) == reverse(circshift(collect(hull2.hull),Int(length(hull2.hull)/2-1))) 
        hullcollinear2 = jarvismarch(boxcoords; collinear=true, by=by)
        @test collect(hullcollinear2.hull) == [[(first(jrange),i) for i in reverse(irange)]..., [(j,first(irange)) for j in jrange[2:end]]..., [(last(jrange),i) for i in irange[2:end]]..., [(j,last(irange)) for j in reverse(jrange)[2:end-1]]...]
        @test collect(jarvismarch(boxcoords; orientation=CW, collinear=true, by=by).hull) == reverse(circshift(collect(hullcollinear2.hull),Int(length(hullcollinear2.hull)/2-1))) 
    end
end