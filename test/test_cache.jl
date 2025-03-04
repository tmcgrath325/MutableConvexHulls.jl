using Random
using PairedLinkedLists: getnode, SkipListCache
using MutableConvexHulls: ChanHullCache, chanhullsidentical

function build_random_hull(H::Type{<:AbstractConvexHull}, n::Int=1024, m::Int=4; subhullcaches=true, orientation=CCW, collinear::Bool=false, sortedby=identity)
    coords = [(rand(), randn()) for i=1:n]

    h = H{eltype(coords)}(orientation, collinear, sortedby)
    h.cache = ChanHullCache{eltype(coords)}()
    if subhullcaches
        h.subhulls[1].points.cache = SkipListCache{eltype(coords)}()
    end
    for i=1:Int(ceil(n/m))
        mergepoints!(h, coords[(i-1)*m+1:min(i*m,n)])
        popidx = rand(1:length(h.hull))
        removepoint!(h, getnode(h.hull, popidx))
    end
    return h
end

@testset "cache" begin
    for H in [ChanConvexHull, ChanLowerConvexHull, ChanUpperConvexHull]
        for o in [CCW, CW]
            for c in [false, true]
                for sb in [identity, x -> (x[1], -x[2])]
                    for shcaches in [false, true]
                        for i=1:10
                            h = build_random_hull(H; subhullcaches=shcaches, orientation=o, collinear=c, sortedby=sb)
                            h2 = MutableConvexHulls.copyfromcache(h)
                            @test h == h2
                            if shcaches
                                @test chanhullsidentical(h, h2)
                            else
                                @test !chanhullsidentical(h, h2)
                            end
                        end
                    end
                end
            end
        end
    end
end