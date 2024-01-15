using Revise
using SpinAdaptedSecondQuantization

abstract type PositronOrbital <: SASQ.Orbital end
abstract type OccPositronOrbital <: PositronOrbital end
abstract type VirPositronOrbital <: PositronOrbital end

SASQ.getnames(::Type{PositronOrbital}) = "PQRSTUVW"
SASQ.getnames(::Type{OccPositronOrbital}) = "IJKLMNO"
SASQ.getnames(::Type{VirPositronOrbital}) = "ABCDEFG"

SASQ.getshortname(::Type{PositronOrbital}) = "GP"
SASQ.getshortname(::Type{OccPositronOrbital}) = "OP"
SASQ.getshortname(::Type{VirPositronOrbital}) = "VP"

SASQ.set_color(PositronOrbital, :purple)
SASQ.set_color(OccPositronOrbital, :blue)
SASQ.set_color(VirPositronOrbital, :red)
