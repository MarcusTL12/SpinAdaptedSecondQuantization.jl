export GeneralOrbital, OccupiedOrbital, VirtualOrbital
export gen, occ, vir, isocc, isvir

"""
    GeneralOrbital

Top level supertype of all orbital spaces.
Any orbital subspace is a subtype of this type.
"""
abstract type GeneralOrbital end

abstract type OccupiedOrbital <: GeneralOrbital end
abstract type VirtualOrbital <: GeneralOrbital end

struct MOIndex{S<:GeneralOrbital}
    name::String
    space::Type{S}
end

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

isocc(::MOIndex{S}) where {S<:GeneralOrbital} =
    (S <: OccupiedOrbital) && !(S <: VirtualOrbital)
isvir(::MOIndex{S}) where {S<:GeneralOrbital} =
    (S <: VirtualOrbital) && !(S <: OccupiedOrbital)
