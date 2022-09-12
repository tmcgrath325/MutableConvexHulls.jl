using PairedLinkedLists: AbstractPairedListNode, AbstractPairedSkipNode, AbstractPairedLinkedList, AbstractPairedSkipList

abstract type AbstractHullNode{T,L,F} <: AbstractPairedListNode{T,L} end
abstract type AbstractHullList{T,F} <: AbstractPairedLinkedList{T} end
abstract type AbstractPointList{T,R,N,F} <: AbstractPairedSkipList{T,F} end


mutable struct PointNode{T,L<:AbstractPointList{T},N<:AbstractHullNode{T}} <: AbstractPairedSkipNode{T,L}
    list::L
    data::T
    prev::PointNode{T,L,N}
    next::PointNode{T,L,N}
    up::PointNode{T,L,N}
    down::PointNode{T,L,N}
    target::Union{N,PointNode{T,L,N}}
    function PointNode{T,L,N}(list::L) where {T,L<:AbstractPointList{T},N<:AbstractHullNode{T}}
        node = new{T,L,N}(list)
        node.next = node
        node.prev = node
        node.target = node
        node.up = node
        node.down = node
        return node
    end
    function PointNode{T,L,N}(list::L, data) where {T,L<:AbstractPointList{T},N<:AbstractHullNode{T}}
        node = new{T,L,N}(list, data)
        node.next = node
        node.prev = node
        node.up = node
        node.down = node
        node.target = node
        return node
    end
end

mutable struct PointList{T,R<:AbstractHullList{T},N<:AbstractHullNode{T},F<:Function} <: AbstractPointList{T,R,N,F}
    len::Int
    nlevels::Int
    skipfactor::Int
    sortedby::F
    target::Union{R,PointList{T,R,N,F}}
    head::PointNode{T,PointList{T,R,N,F},N}
    tail::PointNode{T,PointList{T,R,N,F},N}
    top::PointNode{T, PointList{T,R,N,F},N}
    toptail::PointNode{T, PointList{T,R,N,F},N}
    function PointList{T,R,N,F}(;sortedby::F=identity, skipfactor::Int=2) where {T,R<:AbstractHullList{T},N<:AbstractHullNode{T},F<:Function}
        l = new{T,R,N,F}(0,1,skipfactor,sortedby)
        l.target = l
        l.head = PointNode{T,PointList{T,R,N,F},N}(l)
        l.tail = PointNode{T,PointList{T,R,N,F},N}(l)
        l.top = l.head
        l.toptail = l.tail
        l.sortedby = sortedby
        l.skipfactor = skipfactor
        l.nlevels = 1
        l.head.next = l.tail
        l.tail.prev = l.head
        l.top.next = l.toptail
        l.toptail.prev = l.top
        return l
    end
    function PointList{T,R,N,F}(target::PointList{T}; sortedby::F=identity, skipfactor::Int=2) where {T,R<:AbstractHullList{T},N<:AbstractHullNode{T},F<:Function}
        l = new{T,R,N,F}(0,1,skipfactor,sortedby,target)
        l.head = PointNode{T,PointList{T,R,N,F},N}(l)
        l.tail = PointNode{T,PointList{T,R,N,F},N}(l)
        l.top = l.head
        l.toptail = l.tail
        l.sortedby = sortedby
        l.skipfactor = skipfactor
        l.nlevels = 1
        l.head.next = l.tail
        l.tail.prev = l.head
        l.top.next = l.toptail
        l.toptail.prev = l.top
        return l
    end
end
mutable struct HullNode{T,L<:AbstractHullList{T},F<:Function} <: AbstractHullNode{T,L,F}
    list::L
    data::T
    prev::HullNode{T,L,F}
    next::HullNode{T,L,F}
    target::Union{HullNode{T,L,F},PointNode{T,PointList{T,L,HullNode{T,L,F},F},HullNode{T,L,F}}}
    function HullNode{T,L,F}(list::L) where {T,L<:AbstractHullList{T},F<:Function}
        node = new{T,L,F}(list)
        node.next = node
        node.prev = node
        node.target = node
        return node
    end
    function HullNode{T,L,F}(list::L, data) where {T,L<:AbstractHullList{T},F<:Function}
        node = new{T,L,F}(list, data)
        node.next = node
        node.prev = node
        node.target = node
        return node
    end
end

mutable struct HullList{T,F<:Function} <: AbstractHullList{T,F}
    len::Int
    target::Union{HullList{T,F}, PointList{T,HullList{T,F},HullNode{T,HullList{T,F},F},F}}
    head::HullNode{T,HullList{T,F},F}
    tail::HullNode{T,HullList{T,F},F}
    function HullList{T,F}() where {T,F<:Function}
        l = new{T,F}(0)
        l.target = l
        l.head = HullNode{T,HullList{T,F},F}(l)
        l.tail = HullNode{T,HullList{T,F},F}(l)
        l.head.next = l.tail
        l.tail.prev = l.head
        return l
    end
    function HullList{T,F}(target::PointList{T}) where {T,F<:Function}
        l = new{T,F}(0, target)
        l.head = HullNode{T,HullList{T,F},F}(l)
        l.tail = HullNode{T,HullList{T,F},F}(l)
        l.head.next = l.tail
        l.tail.prev = l.head
        return l
    end
end

HullList{T}(;sortedby::F=identity) where {T,F<:Function} = HullList{T,F}()
function HullList{T,F}(elts...) where {T,F}
    l = HullList{T,F}()
    for elt in elts
        push!(l, elt)
    end
    return l
end

PointList{T}(;sortedby::F=identity, kwargs...) where {T,F<:Function} = PointList{T,HullList{T,F},HullNode{T,HullList{T,F},F},F}(;sortedby=sortedby,kwargs...)
function PointList{T,R,N,F}(elts...) where {T,R,N,F}
    l = PointList{T,R,N,F}()
    for elt in elts
        push!(l, elt)
    end
    return l
end



PairedLinkedLists.nodetype(::Type{HullList{T,F}}) where {T,F} = HullNode{T,HullList{T,F},F}
PairedLinkedLists.nodetype(::Type{PointList{T,R,N,F}}) where {T,R,N,F} = PointNode{T,PointList{T,R,N,F},N}

function PairedLinkedLists.addtarget!(list::L, target::R) where {L<:Union{HullList, PointList},R<:Union{HullList,PointList}}
    if hastarget(list)     # remove existing targets
        removetarget!(list)
    end
    if hastarget(target)
        removetarget!(target)
    end
    list.target = target
    target.target = list
    return list
end

function PairedLinkedLists.addtarget!(node::N, target::R) where {N<:Union{HullNode, PointNode},R<:Union{HullNode,PointNode}}
    node.list.target === target.list || throw(ArgumentError("The provided node must belong to paired list."))
    if hastarget(node)     # remove existing targets
        removetarget!(node)
    end
    if hastarget(target)
        removetarget!(target)
    end
    node.target = target
    target.target = node
    return node
end
