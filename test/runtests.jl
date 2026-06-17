using MutableConvexHulls
using PairedLinkedLists
using Test
using Random
using Logging
using Aqua
using ExplicitImports
const MCH = MutableConvexHulls

Random.seed!(1234)

tests = [
         "orientation",
         "monotonechain",
         "jarvismarch",
         "convexhull",
         "chanhull",
         "cache",
         "testcases",
         "generic_axes",
        ]

@testset "MutableConvexHulls" begin

include(joinpath(dirname(@__FILE__), "test_funs.jl"))

for t in tests
    fp = joinpath(dirname(@__FILE__), "test_$t.jl")
    println("$fp ...")
    include(fp)
end

@testset "Aqua" begin
    Aqua.test_all(MutableConvexHulls)
end

@testset "ExplicitImports" begin
    test_explicit_imports(MutableConvexHulls;
                          all_explicit_imports_are_public   = false,
                          all_qualified_accesses_are_public = false)
end

end

