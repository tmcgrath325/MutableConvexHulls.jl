using MutableConvexHulls
using PairedLinkedLists
using Test

tests = ["orientation",
         "monotonechain",
         "jarvismarch",
         "inside",
        ]

@testset "MutableConvexHulls" begin

for t in tests
    fp = joinpath(dirname(@__FILE__), "test_$t.jl")
    println("$fp ...")
    include(fp)
end

end # @testset

