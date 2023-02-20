export GeneralOrbital, OccupiedOrbital, VirtualOrbital
export general, occupied, virtual, isoccupied, isvirtual, space

# Note: perhaps redefine this as Union{OccupiedOrbital,VirtualOrbital}
"""
    GeneralOrbital

Top level supertype of all orbital spaces.
Any orbital subspace is a subtype of this type.
"""
abstract type GeneralOrbital end

"""
    OccupiedOrbital

Type representing all occupied orbitals.
"""
abstract type OccupiedOrbital <: GeneralOrbital end

"""
    VirtualOrbital

Type representing all virtual orbitals.
"""
abstract type VirtualOrbital <: GeneralOrbital end

# Defining relations to sort indices

# A space is not less than itself
Base.isless(::Type{S}, ::Type{S}) where {S<:GeneralOrbital} = false
Base.isless(::Type{S1}, ::Type{S2}) where
{S1<:GeneralOrbital,S2<:GeneralOrbital} = !(S2 < S1)

# A subspace of a space is considered "greater" than the parent space
# making it come later when sorting
Base.isless(::Type{S1}, ::Type{S2}) where {S2<:GeneralOrbital,S1<:S2} = false
# Base.isless(::Type{S1}, ::Type{S2}) where {S1<:GeneralOrbital,S2<:S1} = true

# Defining occupied orbitals to come before virtuals.
# This is up for debate
Base.isless(::Type{OccupiedOrbital}, ::Type{OccupiedOrbital}) = false
# Base.isless(::Type{VirtualOrbital}, ::Type{VirtualOrbital}) = false

Base.isless(::Type{VirtualOrbital}, ::Type{OccupiedOrbital}) = false
# Base.isless(::Type{OccupiedOrbital}, ::Type{VirtualOrbital}) = true

function is_strict_subspace(::Type{S1}, ::Type{S2}) where
{S1<:GeneralOrbital,S2<:GeneralOrbital}
    S1 != S2 && S1 <: S2
end

function getnames(::Type{S}) where {S<:GeneralOrbital}
    throw("getnames not implemented for space $S")
end

getnames(::Type{GeneralOrbital}) = "pqrstuv"
getnames(::Type{OccupiedOrbital}) = "ijklmno"
getnames(::Type{VirtualOrbital}) = "abcdefg"

getshortname(::Type{GeneralOrbital}) = "G"
getshortname(::Type{OccupiedOrbital}) = "O"
getshortname(::Type{VirtualOrbital}) = "V"

subscript(i) = join(Char(0x2080 + d) for d in reverse!(digits(i)))

function getname(::Type{S}, i::Int) where {S<:GeneralOrbital}
    names = getnames(S)

    name = names[(i-1)%length(names)+1]

    extraind = (i - 1) ÷ length(names)
    if extraind == 0
        name
    else
        name * subscript(extraind)
    end
end

"""
    MOIndex
"""
struct MOIndex
    index::Int
end

function getname(i::MOIndex)
    getname(GeneralOrbital, i.index)
end

general(index) = MOIndex(index)

function Base.show(io::IO, i::MOIndex)
    print(io, getname(i))
end

Base.isdisjoint(::Type{S1}, ::Type{S2}) where
{S1<:GeneralOrbital,S2<:GeneralOrbital} = typeintersect(S1, S2) == Union{}

function Base.isless(p::MOIndex, q::MOIndex)
    p.index < q.index
end

function exchange_index(p::MOIndex, mapping)
    for (old, new) in mapping
        if p == old
            return new
        end
    end
    p
end

# indices must be sorted
function next_free_index(indices)
    i = 1
    for p in indices
        if p.index == i
            i += 1
        end
    end
    MOIndex(i)
end

# utility functions for getting which sets compose each other

# by default spaces can not be added together to form a new space
function add_spaces(::Type{S1}, ::Type{S2}) where
{S1<:GeneralOrbital,S2<:GeneralOrbital}
    nothing
end

function add_spaces(::Type{OccupiedOrbital}, ::Type{VirtualOrbital})
    GeneralOrbital
end

function add_spaces(::Type{VirtualOrbital}, ::Type{OccupiedOrbital})
    GeneralOrbital
end

# By default taking the difference of spaces does not form a new space
function diff_spaces(::Type{S1}, ::Type{S2}) where
{S1<:GeneralOrbital,S2<:GeneralOrbital}
    nothing
end

function diff_spaces(::Type{GeneralOrbital}, ::Type{OccupiedOrbital})
    VirtualOrbital
end

function diff_spaces(::Type{GeneralOrbital}, ::Type{VirtualOrbital})
    OccupiedOrbital
end
