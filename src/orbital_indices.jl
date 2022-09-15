export GeneralOrbital, OccupiedOrbital, VirtualOrbital
export gen, occ, vir, isocc, isvir, space

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

# A subspace of a space is considered "greater" than the parent space
# making it come later when sorting
Base.isless(::Type{S1}, ::Type{S2}) where {S2<:GeneralOrbital,S1<:S2} = false
Base.isless(::Type{S1}, ::Type{S2}) where {S1<:GeneralOrbital,S2<:S1} = true

# Defining occupied orbitals to come before virtuals.
# This is up for debate
Base.isless(::Type{OccupiedOrbital}, ::Type{OccupiedOrbital}) = false
Base.isless(::Type{VirtualOrbital}, ::Type{VirtualOrbital}) = false

Base.isless(::Type{VirtualOrbital}, ::Type{OccupiedOrbital}) = false
Base.isless(::Type{OccupiedOrbital}, ::Type{VirtualOrbital}) = true


"""
    MOIndex
"""
struct MOIndex{S<:GeneralOrbital}
    name::String
    function MOIndex(name, ::Type{S}) where {S<:GeneralOrbital}
        new{S}(name)
    end
end

space(::MOIndex{S}) where {S<:GeneralOrbital} = S

gen(name) = MOIndex(name, GeneralOrbital)
occ(name) = MOIndex(name, OccupiedOrbital)
vir(name) = MOIndex(name, VirtualOrbital)

function Base.show(io::IO, i::MOIndex{GeneralOrbital})
    print(io, i.name)
end

function Base.show(io::IO, i::MOIndex{OccupiedOrbital})
    print(io, "\x1b[92m", i.name, "\x1b[39m")
end

function Base.show(io::IO, i::MOIndex{VirtualOrbital})
    print(io, "\x1b[36m", i.name, "\x1b[39m")
end

Base.isdisjoint(::MOIndex{S1}, ::MOIndex{S2}) where
{S1<:GeneralOrbital,S2<:GeneralOrbital} = typeintersect(S1, S2) == Union{}

isocc(::MOIndex{S}) where {S<:GeneralOrbital} = S <: OccupiedOrbital
isvir(::MOIndex{S}) where {S<:GeneralOrbital} = S <: VirtualOrbital

function Base.isless(p::MOIndex{S1}, q::MOIndex{S2}) where
{S1<:GeneralOrbital,S2<:GeneralOrbital}
    (S1, p.name) < (S2, q.name)
end
