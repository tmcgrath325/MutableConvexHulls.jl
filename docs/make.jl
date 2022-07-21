using MutableConvexHulls
using Documenter

DocMeta.setdocmeta!(MutableConvexHulls, :DocTestSetup, :(using MutableConvexHulls); recursive=true)

makedocs(;
    modules=[MutableConvexHulls],
    authors="Tom McGrath <tmcgrath325@gmail.com> and contributors",
    repo="https://github.com/tmcgrath325/MutableConvexHulls.jl/blob/{commit}{path}#{line}",
    sitename="MutableConvexHulls.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://tmcgrath325.github.io/MutableConvexHulls.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/tmcgrath325/MutableConvexHulls.jl",
    devbranch="main",
)
