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

function getnames(::Type{S}) where {S<:GeneralOrbital}
    throw("getnames not implemented for space $S")
end

getnames(::Type{GeneralOrbital}) = "pqrstuv"
getnames(::Type{OccupiedOrbital}) = "ijklmno"
getnames(::Type{VirtualOrbital}) = "abcdefg"

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
struct MOIndex{S<:GeneralOrbital}
    index::Int
    function MOIndex(index, ::Type{S}) where {S<:GeneralOrbital}
        new{S}(index)
    end
end

space(::MOIndex{S}) where {S<:GeneralOrbital} = S

function getname(i::MOIndex{S}) where {S<:GeneralOrbital}
    getname(S, i.index)
end

general(index) = MOIndex(index, GeneralOrbital)
occupied(index) = MOIndex(index, OccupiedOrbital)
virtual(index) = MOIndex(index, VirtualOrbital)

function Base.show(io::IO, i::MOIndex{GeneralOrbital})
    print(io, getname(i))
end

function Base.show(io::IO, i::MOIndex{OccupiedOrbital})
    print(io, "\x1b[92m", getname(i), "\x1b[39m")
end

function Base.show(io::IO, i::MOIndex{VirtualOrbital})
    print(io, "\x1b[36m", getname(i), "\x1b[39m")
end

Base.isdisjoint(::MOIndex{S1}, ::MOIndex{S2}) where
{S1<:GeneralOrbital,S2<:GeneralOrbital} = typeintersect(S1, S2) == Union{}

isoccupied(::MOIndex{S}) where {S<:GeneralOrbital} = S <: OccupiedOrbital
isvirtual(::MOIndex{S}) where {S<:GeneralOrbital} = S <: VirtualOrbital

function Base.isless(p::MOIndex{S1}, q::MOIndex{S2}) where
{S1<:GeneralOrbital,S2<:GeneralOrbital}
    (S1, p.index) < (S2, q.index)
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
function next_free_index(indices, ::Type{S}) where {S<:GeneralOrbital}
    i = 1
    for p in indices
        if space(p) == S
            if p.index == i
                i += 1
            end
        end
    end
    MOIndex(i, S)
end

function next_free_index(indices, ::MOIndex{S}) where {S<:GeneralOrbital}
    next_free_index(indices, S)
end
