@testset "convex hulls, unique points" begin
    boxcoords = [(i,j) for i in 1:10 for j in 1:10]
    by = x -> (x[1], -x[2])
    n = 10
    
    hulltestset("lower convex hull", n, by, boxcoords, MutableLowerConvexHull, lower_jarvismarch)
    hulltestset("upper convex hull", n, by, boxcoords, MutableUpperConvexHull, upper_jarvismarch)
    hulltestset("convex hull",       n, by, boxcoords, MutableConvexHull,      jarvismarch)
end

@testset "convex hulls, duplicate points" begin
    boxcoords = [(i,j) for i in 1:10 for j in 1:10]
    boxcoords = [boxcoords..., boxcoords..., boxcoords...]
    by = x -> (x[1], -x[2])
    n = 10
    
    hulltestset("lower convex hull", n, by, boxcoords, MutableLowerConvexHull, lower_jarvismarch)
    hulltestset("upper convex hull", n, by, boxcoords, MutableUpperConvexHull, upper_jarvismarch)
    hulltestset("convex hull",       n, by, boxcoords, MutableConvexHull,      jarvismarch)
end

@testset "convex hulls, random data with duplicates" begin
    coords = [(randn(),randn()) for i in 1:10 for j in 1:10]
    by = x -> (x[1], -x[2])
    n = 10
    
    hulltestset("lower convex hull", n, by, coords, MutableLowerConvexHull, lower_jarvismarch)
    hulltestset("upper convex hull", n, by, coords, MutableUpperConvexHull, upper_jarvismarch)
    hulltestset("convex hull",       n, by, coords, MutableConvexHull,      jarvismarch)
end