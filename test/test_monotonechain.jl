# using MutableConvexHulls: monotonechain!, lower_monotonechain!, upper_monotonechain!

@testset "Monotone Chain Algorithm" begin
    irange = 1:10
    jrange = 1:10
    boxcoords = [(i,j) for i in irange for j in jrange]
    dupcoords = [boxcoords..., boxcoords...]
    by = x -> (x[1], -x[2])

    @testset "Lower Monotone Chain" begin
        # standard-sorted coords
        lower = lower_monotonechain(boxcoords) # ; orientation = CCW, collinear = false
        @test lower == lower_monotonechain(dupcoords)
        @test collect(lower.hull) == collect(lower_monotonechain(boxcoords).hull) == [first(boxcoords), (last(boxcoords)[1], first(boxcoords)[2]), last(boxcoords)]
        @test collect(lower_monotonechain(boxcoords; orientation=CW).hull) == reverse(collect(lower.hull))
        lowercollinear = lower_monotonechain(boxcoords; collinear = true)
        @test lowercollinear == lower_monotonechain(dupcoords; collinear = true)
        @test collect(lowercollinear.hull) == collect(lower_monotonechain(boxcoords; collinear=true).hull) == [[(i,first(jrange)) for i in irange]..., [(last(irange), j) for j in jrange[2:end]]...]
        @test collect(lower_monotonechain(boxcoords; orientation=CW, collinear=true).hull) == reverse(collect(lowercollinear.hull))

        # alternate-sorting coords
        lower2 = lower_monotonechain(boxcoords; sortedby=by) # ; orientation = CCW, collinear = false
        @test lower2 == lower_monotonechain(dupcoords; sortedby=by)
        @test collect(lower2.hull) == [(first(boxcoords)[1], last(boxcoords)[2]), first(boxcoords), (last(boxcoords)[1], first(boxcoords)[2])]
        @test collect(lower_monotonechain(boxcoords; orientation=CW, sortedby=by).hull) == reverse(collect(lower2.hull))
        lowercollinear2 = lower_monotonechain(boxcoords; collinear=true, sortedby=by)
        @test lowercollinear2 == lower_monotonechain(dupcoords; collinear=true, sortedby=by)
        @test collect(lowercollinear2.hull) == [[(first(jrange),i) for i in reverse(irange)]..., [(j,first(irange)) for j in jrange[2:end]]...]
        @test collect(lower_monotonechain(boxcoords; orientation=CW, collinear=true, sortedby=by).hull) == reverse(collect(lowercollinear2.hull))
    end

    @testset "Upper Monotone Chain" begin
        # standard-sorted coords
        upper = upper_monotonechain(boxcoords) # ; orientation = CCW, collinear = false
        @test upper == upper_monotonechain(dupcoords)
        @test collect(upper.hull) == collect(upper_monotonechain(boxcoords).hull) == [last(boxcoords), (first(boxcoords)[1], last(boxcoords)[2]), first(boxcoords)]
        @test collect(upper_monotonechain(boxcoords; orientation=CW).hull) == reverse(collect(upper.hull))
        uppercollinear = upper_monotonechain(boxcoords; collinear = true)
        @test uppercollinear == upper_monotonechain(dupcoords; collinear = true)
        @test collect(uppercollinear.hull) == collect(upper_monotonechain(boxcoords; collinear=true).hull) == [[(i,last(jrange)) for i in reverse(irange)]..., [(first(irange), j) for j in reverse(jrange)[2:end]]...]
        @test collect(upper_monotonechain(boxcoords; orientation=CW, collinear=true).hull) == reverse(collect(uppercollinear.hull))

        # alternate-sorting coords
        upper2 = upper_monotonechain(boxcoords; sortedby=by) # ; orientation = CCW, collinear = false
        @test upper2 == upper_monotonechain(dupcoords; sortedby=by)
        @test collect(upper2.hull) == [(last(boxcoords)[1], first(boxcoords)[2]), last(boxcoords), (first(boxcoords)[1], last(boxcoords)[2])]
        @test collect(upper_monotonechain(boxcoords; orientation=CW, sortedby=by).hull) == reverse(collect(upper2.hull))
        uppercollinear2 = upper_monotonechain(boxcoords; collinear=true, sortedby=by)
        @test uppercollinear2 == upper_monotonechain(dupcoords; collinear=true, sortedby=by)
        @test collect(uppercollinear2.hull) == [[(last(jrange),i) for i in irange]..., [(j,last(irange)) for j in reverse(jrange)[2:end]]...]
        @test collect(upper_monotonechain(boxcoords; orientation=CW, collinear=true, sortedby=by).hull) == reverse(collect(uppercollinear2.hull))
    end

    @testset "Full Monotone Chain" begin
        # standard-sorted coords
        hull = monotonechain(boxcoords) # ; orientation = CCW, collinear = false
        @test hull == monotonechain(dupcoords)
        @test collect(hull.hull) == collect(monotonechain(boxcoords).hull) == [first(boxcoords), (last(boxcoords)[1], first(boxcoords)[2]), last(boxcoords), (first(boxcoords)[1], last(boxcoords)[2])]
        @test collect(monotonechain(boxcoords; orientation=CW).hull) == reverse(circshift(collect(hull.hull),Int(length(hull.hull)/2-1))) 
        hullcollinear = monotonechain(boxcoords; collinear = true)
        @test hullcollinear == monotonechain(dupcoords; collinear = true)
        @test collect(hullcollinear.hull) == collect(monotonechain(boxcoords; collinear=true).hull) == [[(i,first(jrange)) for i in irange]..., [(last(irange), j) for j in jrange[2:end]]..., [(i,last(jrange)) for i in reverse(irange)[2:end]]..., [(first(irange), j) for j in reverse(jrange)[2:end-1]]...]
        @test collect(monotonechain(boxcoords; orientation=CW, collinear=true).hull) == reverse(circshift(collect(hullcollinear.hull), Int(length(hullcollinear.hull)/2)-1))

        # alternate-sorting coords
        hull2 = monotonechain(boxcoords; sortedby=by) # ; orientation = CCW, collinear = false
        @test hull2 == monotonechain(dupcoords; sortedby=by)
        @test collect(hull2.hull) == [(first(boxcoords)[1], last(boxcoords)[2]), first(boxcoords), (last(boxcoords)[1], first(boxcoords)[2]), last(boxcoords)]
        @test collect(monotonechain(boxcoords; orientation=CW, sortedby=by).hull) == reverse(circshift(collect(hull2.hull),Int(length(hull2.hull)/2-1))) 
        hullcollinear2 = monotonechain(boxcoords; collinear=true, sortedby=by)
        @test hullcollinear2 == monotonechain(dupcoords; collinear=true, sortedby=by)
        @test collect(hullcollinear2.hull) == [[(first(jrange),i) for i in reverse(irange)]..., [(j,first(irange)) for j in jrange[2:end]]..., [(last(jrange),i) for i in irange[2:end]]..., [(j,last(irange)) for j in reverse(jrange)[2:end-1]]...]
        @test collect(monotonechain(boxcoords; orientation=CW, collinear=true, sortedby=by).hull) == reverse(circshift(collect(hullcollinear2.hull),Int(length(hullcollinear2.hull)/2-1))) 
    end
end