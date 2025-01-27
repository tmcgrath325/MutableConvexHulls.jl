@testset "chan hulls, unique points" begin
    boxcoords = [(i,j) for i in 1:10 for j in 1:10]
    by = x -> (x[1], -x[2])
    n = 10
    
    chanhulltestset("lower chan hull", n, by, boxcoords, ChanLowerConvexHull, lower_jarvismarch)
    chanhulltestset("upper chan hull", n, by, boxcoords, ChanUpperConvexHull, upper_jarvismarch)
    chanhulltestset("chan hull",       n, by, boxcoords, ChanConvexHull,      jarvismarch)
end

@testset "chan hulls, duplicate points" begin
    boxcoords = [(i,j) for i in 1:10 for j in 1:10]
    boxcoords = [boxcoords..., boxcoords..., boxcoords...]
    by = x -> (x[1], -x[2])
    n = 10
    
    chanhulltestset("lower chan hull", n, by, boxcoords, ChanLowerConvexHull, lower_jarvismarch)
    chanhulltestset("upper chan hull", n, by, boxcoords, ChanUpperConvexHull, upper_jarvismarch)
    chanhulltestset("chan hull",       n, by, boxcoords, ChanConvexHull,      jarvismarch)
end

@testset "chan hulls, random data with duplicates" begin
    coords = [(randn(),randn()) for i in 1:10 for j in 1:10]
    coords = [coords..., coords..., coords...]
    by = x -> (x[1], -x[2])
    n = 10
    
    chanhulltestset("lower chan hull", n, by, coords, ChanLowerConvexHull, lower_jarvismarch)
    chanhulltestset("upper chan hull", n, by, coords, ChanUpperConvexHull, upper_jarvismarch)
    chanhulltestset("chan hull",       n, by, coords, ChanConvexHull,      jarvismarch)
end