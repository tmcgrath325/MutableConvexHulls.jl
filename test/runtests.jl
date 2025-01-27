using MutableConvexHulls
using PairedLinkedLists
using Test
using Random
const MCH = MutableConvexHulls

Random.seed!(1234)

tests = ["orientation",
         "monotonechain",
         "jarvismarch",
         "convexhull",
         "chanhull",
         "testcases"
        ]

@testset "MutableConvexHulls" begin

include(joinpath(dirname(@__FILE__), "test_funs.jl"))

for t in tests
    fp = joinpath(dirname(@__FILE__), "test_$t.jl")
    println("$fp ...")
    include(fp)
end

end

