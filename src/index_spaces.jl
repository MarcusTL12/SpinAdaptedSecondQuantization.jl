const index_ids::Dict{Any,Int} = Dict()

const index_shortnames::Vector{String} = String[]

const index_names::Vector{String} = String[]

function new_space(id, shortname::String, names::String)
    new_ind = length(index_ids) + 1
    index_ids[id] = new_ind

    push!(index_shortnames, shortname)
    push!(index_names, names)

    new_ind
end

const GeneralIndex = new_space(:GeneralIndex, "GI", "")

function getshortname(s::Int)
    index_shortnames[s]
end

const space_sums::Dict{NTuple{2,Int},Int} = Dict()
const space_diff::Dict{NTuple{2,Int},Int} = Dict()

function add_space_sum(a::Int, b::Int, c::Int)
    a, b = a < b ? (a, b) : (b, a)
    space_sums[(a, b)] = c
    space_diff[(c, a)] = b
    space_diff[(c, b)] = a
end

const space_intersections::Dict{NTuple{2,Int},Int} = Dict()

function add_space_intersection(a::Int, b::Int, c::Int)
    a, b = a < b ? (a, b) : (b, a)
    space_intersections[(a, b)] = c
end

function Base.intersect(a::Int, b::Int)
    if a == b
        a
    else
        get(space_intersections, a < b ? (a, b) : (b, a), nothing)
    end
end

function Base.isdisjoint(a::Int, b::Int)
    if a == b
        return false
    end

    a, b = a < b ? (a, b) : (b, a)
    !haskey(space_intersections, (a, b))
end

function Base.:⊆(a::Int, b::Int)
    intersect(a, b) == a
end

function Base.:⊊(a::Int, b::Int)
    a != b && a ⊆ b
end

function add_spaces(a::Int, b::Int)
    a, b = a < b ? (a, b) : (b, a)
    get(space_sums, (a, b), nothing)
end

function diff_spaces(a::Int, b::Int)
    get(space_diff, (a, b), nothing)
end

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

colors::Dict{Int,Union{Int,Int}} = Dict()

function getname(io::IO, s::Int, i::Int)
    if s == GeneralIndex
        subscript(io, i)
        return
    end

    names = index_names[s]

    print(io, names[(i-1)%length(names)+1])

    extraind = (i - 1) ÷ length(names)
    if extraind != 0
        subscript(io, extraind)
    end
end

"""
    Constraints = SortedDict{Int,Int}

Type alias for container of Index constraints
"""
const Constraints = SortedDict{Int,Int}

# Lets you call the constraints as a function instead of explicitly calling get
function (constraints::Constraints)(p::Int)
    get(constraints, p, GeneralIndex)
end

const IndexTranslation = Dict{Int,Tuple{Int,Int}}

function (translation::IndexTranslation)(p::Int)
    get(translation, p, (:GeneralIndex, p))
end

export translate
function translate(translations...)
    counts = Dict{Int,Int}()
    translation = IndexTranslation()

    for (S, inds) in translations
        for p in inds
            counts[S] = get(counts, S, 0) + 1
            translation[p] = (S, counts[S])
        end
    end

    translation
end

function getname(io::IO, constraints::Constraints,
    translation::IndexTranslation, i::Int)
    do_color = index_color && haskey(colors, constraints(i))

    if do_color
        print(io, Base.text_colors[colors[constraints(i)]])
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
do_index_translation::Bool = true

export enable_color, disable_color, set_color,
    enable_index_translation, disable_index_translation

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
    set_color(S, color=default_color(S))

Sets the color of a given orbital space to the specified color, or the
default if none is given. Allowed colors are given by the
[Base.text_colors](https://github.com/JuliaLang/julia/blob/\
17cfb8e65ead377bf1b4598d8a9869144142c84e/base/util.jl#L5-L34)
dict. This includes integers in the range 0-255 which the corresponding colors
can be found in this table on
[Wikipedia](https://en.wikipedia.org/wiki/ANSI_escape_code#8-bit).
"""
function set_color(s::Int, color)
    colors[s] = color
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
function disable_color(s::Int)
    delete!(colors, s)
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
function enable_index_translation()
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
