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

getnames(::Type{GeneralOrbital}) = "pqrstuvw"
getnames(::Type{OccupiedOrbital}) = "ijklmno"
getnames(::Type{VirtualOrbital}) = "abcdefg"

getshortname(::Type{GeneralOrbital}) = "G"
getshortname(::Type{OccupiedOrbital}) = "O"
getshortname(::Type{VirtualOrbital}) = "V"

function subscript(io::IO, i)
    for d in reverse!(digits(i))
        print(io, Char(0x2080 + d))
    end
end

function subscript(i)
    io = IOBuffer()
    subscript(io, i)
    String(take!(io))
end

function getname(io::IO, ::Type{S}, i::Int) where {S<:GeneralOrbital}
    names = getnames(S)

    print(io, names[(i-1)%length(names)+1])

    extraind = (i - 1) ÷ length(names)
    if extraind != 0
        subscript(io, extraind)
    end
end

function getname(::Type{S}, i::Int) where {S<:GeneralOrbital}
    io = IOBuffer()
    getname(io, S, i)
    String(take!(io))
end

function print_mo_index(io::IO, p)
    getname(io, GeneralOrbital, p)
end

function print_mo_index(io::IO, indices...)
    for p in indices
        print_mo_index(io, p)
    end
end

function print_mo_index(p)
    getname(GeneralOrbital, p)
end

Base.isdisjoint(::Type{S1}, ::Type{S2}) where
{S1<:GeneralOrbital,S2<:GeneralOrbital} = typeintersect(S1, S2) == Union{}

function exchange_index(p::Int, mapping)
    for (old, new) in mapping
        if p == old
            return new
        end
    end
    p
end

# indices must be sorted
function next_free_index(indices)
    @assert issorted(indices)
    i = 1
    for p in indices
        if i < p
            return p
        end
        i += 1
    end
    i
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
