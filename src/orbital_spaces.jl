export GeneralOrbital, OccupiedOrbital, VirtualOrbital
export general, isoccupied, isvirtual

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

default_color(::Type{GeneralOrbital}) = :nothing
default_color(::Type{OccupiedOrbital}) = :light_green
default_color(::Type{VirtualOrbital}) = :cyan

colors::Dict{Type,Union{Symbol,Int}} = Dict{Type,Union{Symbol,Int}}([
    GeneralOrbital => default_color(GeneralOrbital),
    OccupiedOrbital => default_color(OccupiedOrbital),
    VirtualOrbital => default_color(VirtualOrbital),
])

function getname(io::IO, ::Type{S}, i::Int) where {S<:GeneralOrbital}
    names = getnames(S)

    print(io, names[(i-1)%length(names)+1])

    extraind = (i - 1) ÷ length(names)
    if extraind != 0
        subscript(io, extraind)
    end
end

const IndexTranslation = Dict{Int,Tuple{Type,Int}}

function (translation::IndexTranslation)(p::Int)
    get(translation, p, (GeneralOrbital, p))
end

function getname(io::IO, constraints::Constraints,
    translation::IndexTranslation, i::Int)
    do_color = index_color &&
               (is_strict_subspace(constraints(i), translation(i)[1]) ||
                color_translated)

    if do_color
        print(io, Base.text_colors[get(colors, constraints(i), :nothing)])
    end

    getname(io, translation(i)...)

    if do_color
        print(io, "\x1b[39m")
    end
end

function print_mo_index(io::IO, constraints::Constraints,
    translation::IndexTranslation, p)
    getname(io, constraints, translation, p)
end

function print_mo_index(io::IO, constraints::Constraints,
    translation::IndexTranslation, indices...)
    for p in indices
        print_mo_index(io, constraints, translation, p)
    end
end

index_color::Bool = true
color_translated::Bool = false
do_index_translation::Bool = true

export enable_color, disable_color, set_color,
    enable_color_translated, disable_color_translated,
    enable_index_tranlation, disable_index_translation

"""
    enable_color()

Enables the coloring of MO-indices to indicate orbital space constraints.
Color is enabled by default.
The default coloring scheme is given by:
- GeneralOrbital: :nothing
- OccupiedOrbital: :light_green
- VirtualOrbital: :cyan
"""
function enable_color()
    global index_color = true
    nothing
end

"""
    set_color(S, [color=default_color(S)])

Sets the color of a given orbital space to the specified color, or the
default if none is given. Allowed colors are given by the
[Base.text_colors](https://github.com/JuliaLang/julia/blob/\
17cfb8e65ead377bf1b4598d8a9869144142c84e/base/util.jl#L5-L34)
dict. This includes integers in the range 0-255 which the corresponding colors
can be found in this table on
[Wikipedia](https://en.wikipedia.org/wiki/ANSI_escape_code#8-bit).
"""
function set_color(::Type{S}, color=default_color(S)) where
{S<:GeneralOrbital}
    colors[S] = color
    nothing
end

"""
    disable_color()

Disables the coloring of MO-indices to indicate orbital space constraints
globally. This will make all constraints be printed explicitly which can make
some terms a bit long to read. Color is enabled by default.
"""
function disable_color()
    global index_color = false
    nothing
end

"""
    disable_color(S)

Disables the coloring for a specific orbital space. This can be useful
when dealing with a medium-large number of orbital spaces, and having
distinguishable colors is infeasible.
"""
function disable_color(::Type{S}) where {S<:GeneralOrbital}
    delete!(colors, S)
    nothing
end

"""
    enable_index_translation()

Enables the translation of summation indices over occupied and virtual to be
printed as ijkl... and abcd... respectively instead of pqrs...

By default translated indices lose their subspace coloring to reduce redundant
information (unless the index is in an even stricter subspace which requires
coloring nevertheless). To re-enable the coloring see
[`enable_color_translated`](@ref).

This is enabled by default.
"""
function enable_index_tranlation()
    global do_index_translation = true
    nothing
end

"""
    disables_index_translation()

Disables the translation of summation indices over occupied and virtual to be
printed as ijkl... and abcd... respectively instead of pqrs...

See [`enable_index_tranlation`](@ref)
"""
function disable_index_translation()
    global do_index_translation = false
    nothing
end

"""
    enable_color_translated()

Enables the coloring of indices that has been translated.
"""
function enable_color_translated()
    global color_translated = true
    nothing
end

"""
    disable_color_translated()

Disables the coloring of indices that has been translated.
"""
function disable_color_translated()
    global color_translated = false
    nothing
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

function next_free_index(indices)
    sort!(indices)
    i = 1
    for p in indices
        if i == p
            i += 1
        end
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
