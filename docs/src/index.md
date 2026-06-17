```@meta
CurrentModule = MutableConvexHulls
DocTestSetup = :(using MutableConvexHulls)
```

# MutableConvexHulls

MutableConvexHulls provides types and operations for computing and incrementally updating
convex hulls of 2-D point sets. It is designed for use cases where points are added or
removed one at a time and the hull must remain current after each change.

## Installation

MutableConvexHulls is a registered Julia package. Install it with:

```julia
using Pkg
Pkg.add("MutableConvexHulls")
```

Or in the package manager REPL (`]`):

```
pkg> add MutableConvexHulls
```

## Quick start

Build a hull from a collection of points, then add, query, and remove points incrementally:

```jldoctest quickstart
julia> h = monotonechain([(0.0, 0.0), (1.0, 0.0), (0.5, 0.5), (1.0, 1.0), (0.0, 1.0)])
MutableConvexHull{Tuple{Float64, Float64}, typeof(identity)}((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0))

julia> h, expanded = addpoint!(h, (2.0, 0.5));

julia> expanded
true

julia> collect(h)
5-element Vector{Tuple{Float64, Float64}}:
 (0.0, 0.0)
 (1.0, 0.0)
 (2.0, 0.5)
 (1.0, 1.0)
 (0.0, 1.0)

julia> (2.0, 0.5) in h  # on the boundary
true

julia> (3.0, 0.0) in h  # outside
false

julia> h, removed = removepoint!(h, (0.0, 0.0));

julia> removed
true

julia> collect(h)
4-element Vector{Tuple{Float64, Float64}}:
 (0.0, 1.0)
 (1.0, 0.0)
 (2.0, 0.5)
 (1.0, 1.0)
```

## Concepts

### Full, lower, and upper hulls

Three hull variants are available for each hull family:

- **Full hull** ([`MutableConvexHull`](@ref), [`ChanConvexHull`](@ref)) — the complete convex
  hull boundary, listed as a closed polygon.
- **Lower hull** ([`MutableLowerConvexHull`](@ref), [`ChanLowerConvexHull`](@ref)) — the chain
  of vertices along the bottom boundary, from the leftmost to the rightmost point.
- **Upper hull** ([`MutableUpperConvexHull`](@ref), [`ChanUpperConvexHull`](@ref)) — the chain
  of vertices along the top boundary, from the rightmost to the leftmost point.

```jldoctest partialhulls
julia> pts = [(0.0, 0.0), (1.0, -0.5), (2.0, 0.0), (1.5, 1.0), (0.5, 1.0)];

julia> lower_monotonechain(pts)
MutableLowerConvexHull{Tuple{Float64, Float64}, typeof(identity)}((0.0, 0.0), (1.0, -0.5), (2.0, 0.0))

julia> upper_monotonechain(pts)
MutableUpperConvexHull{Tuple{Float64, Float64}, typeof(identity)}((2.0, 0.0), (1.5, 1.0), (0.5, 1.0), (0.0, 0.0))
```

### `MutableConvexHull` vs `ChanConvexHull`

Both families maintain the same hull and support the same operations, but they use different
internal representations:

- **`MutableConvexHull`** stores all points in a single sorted linked list. Each `addpoint!`
  or `removepoint!` call checks whether the new point changes the boundary and, if so,
  reruns the monotone-chain algorithm over the affected segment. This is efficient for hulls
  with frequent point removals or small numbers of points.

- **`ChanConvexHull`** distributes points across a vector of `MutableConvexHull` sub-hulls
  and derives the outer hull by merging their boundaries. Points are never removed from
  sub-hulls during a removal — instead the sub-hull is recomputed — which makes large
  batch additions cheaper. Use `ChanConvexHull` when the point set is large and additions
  dominate.

### Orientation and collinear points

Both families accept two keyword arguments that control how the hull boundary is reported:

- **`orientation`** — `CCW` (default) lists vertices counterclockwise; `CW` lists them
  clockwise.
- **`collinear`** — when `false` (default), collinear points on an edge are not included
  as hull vertices; when `true`, they are.

### Constructors: `monotonechain` vs `jarvismarch`

Two batch constructors build a hull from an existing point collection:

- [`monotonechain`](@ref) uses the [monotone chain algorithm](https://doi.org/10.1016/0020-0190(79)90072-3)
  and is the recommended default.
- [`jarvismarch`](@ref) uses the [gift-wrapping (Jarvis march) algorithm](https://en.wikipedia.org/wiki/Gift_wrapping_algorithm).

Both return the same type with the same interface; subsequent incremental operations use
monotone chain internally regardless of which constructor was used.

### Merging hulls

[`mergehulls!`](@ref) merges one or more hulls into another in place; [`mergehulls`](@ref)
returns a new hull without mutating its arguments:

```jldoctest merge
julia> h1 = monotonechain([(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]);

julia> h2 = monotonechain([(2.0, 0.0), (2.0, 1.0), (1.0, 1.0)]);

julia> mergehulls(h1, h2)
MutableConvexHull{Tuple{Float64, Float64}, typeof(identity)}((0.0, 0.0), (2.0, 0.0), (2.0, 1.0), (0.0, 1.0))
```

### Accessing nodes

Iterating a hull directly yields its vertices as values. For finer-grained access to the
underlying linked-list nodes — for example, to pass a specific node to
[`removepoint!`](@ref) — use the semi-public iterators
[`MutableConvexHulls.HullNodeIterator`](@ref) (over [`HullNode`](@ref)s in `h.hull`) and
[`MutableConvexHulls.PointNodeIterator`](@ref) (over [`PointNode`](@ref)s in `h.points`,
including interior points).

See the [API Reference](@ref) for the full list of exported functions and types.
