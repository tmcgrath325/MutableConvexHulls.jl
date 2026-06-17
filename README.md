# MutableConvexHulls

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://tmcgrath325.github.io/MutableConvexHulls.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://tmcgrath325.github.io/MutableConvexHulls.jl/dev/)
[![Build Status](https://github.com/tmcgrath325/MutableConvexHulls.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/tmcgrath325/MutableConvexHulls.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/tmcgrath325/MutableConvexHulls.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/tmcgrath325/MutableConvexHulls.jl)
[![Aqua QA](https://juliatesting.github.io/Aqua.jl/dev/assets/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

MutableConvexHulls provides a few types that facilitate the calculation of convex hulls from sets of 2D points, as well as supporting efficient updates to the hull upon addition or removal of points.

## Quick start

```julia
using MutableConvexHulls

points = [(0.0,0.0), (1.0,0.0), (1.0,1.0), (0.0,1.0), (0.5,0.5)]
h = monotonechain(points)
collect(h.hull)   # [(0.0,0.0), (1.0,0.0), (1.0,1.0), (0.0,1.0)]

# Add a point — hull updates in place
addpoint!(h, (2.0, 0.5))
collect(h.hull)   # [(0.0,0.0), (1.0,0.0), (2.0,0.5), (1.0,1.0), (0.0,1.0)]

# Membership test (boundary counts as inside when collinear=true)
(2.0, 0.5) in h   # true
(0.5, 0.5) in h   # true  (interior point)

# Remove a point by value — hull updates in place
removepoint!(h, (0.0, 0.0))
collect(h.hull)   # [(0.0,1.0), (1.0,0.0), (2.0,0.5), (1.0,1.0)]
```

`jarvismarch` is an alternative constructor with the same interface. Lower and upper hull variants (`lower_monotonechain`, `upper_jarvismarch`, …) are also available.
