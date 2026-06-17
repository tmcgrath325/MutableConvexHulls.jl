```@meta
CurrentModule = MutableConvexHulls
```

# API Reference

## Hull types

### Abstract types

```@docs
AbstractConvexHull
AbstractChanConvexHull
```

### Full hull

```@docs
MutableConvexHull
ChanConvexHull
```

### Lower hull

```@docs
MutableLowerConvexHull
ChanLowerConvexHull
```

### Upper hull

```@docs
MutableUpperConvexHull
ChanUpperConvexHull
```

## Batch constructors

### Monotone chain

```@docs
monotonechain
lower_monotonechain
upper_monotonechain
```

### Jarvis march

```@docs
jarvismarch
lower_jarvismarch
upper_jarvismarch
```

## Incremental point operations

```@docs
addpoint!
mergepoints!
removepoint!
```

## Merging hulls

```@docs
mergehulls!
mergehulls
```

## Membership testing

```@docs
insidehull
```

## Orientation

```@docs
CCW
CW
```

## Node access

Most users iterate a hull directly to obtain vertex values. When you need to work with the
underlying linked-list nodes — for example, to pass a node to [`removepoint!`](@ref) — use
the iterators and types below.

```@docs
MutableConvexHulls.HullNodeIterator
MutableConvexHulls.PointNodeIterator
HullList
PointList
HullNode
PointNode
```
