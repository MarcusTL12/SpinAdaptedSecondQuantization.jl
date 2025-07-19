const index_ids::Dict{Any,Int} = Dict()

const index_shortnames::Vector{String} = String[]

const index_names::Vector{String} = String[]

function new_space(id, shortname::String, names::String)
    if haskey(index_ids, id)
        return index_ids[id]
    end

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
    add_subspace_relation(c, a)
    add_subspace_relation(c, b)
end

const space_intersections::Dict{NTuple{2,Int},Int} = Dict()

function add_space_intersection(a::Int, b::Int, c::Int)
    a, b = a < b ? (a, b) : (b, a)
    space_intersections[(a, b)] = c

    complete_subspace_tree()
end

function add_subspace_relation(a::Int, c::Int)
    add_space_intersection(a, c, c)
end

function complete_subspace_tree()
    for (_, from) in index_ids
        queue = [from]
        seen = [from]

        while !isempty(queue)
            curnode = popfirst!(queue)
            for (_, s) in index_ids
                if s ⊊ curnode && s ∉ seen
                    if !(s ⊊ from)
                        add_subspace_relation(from, s)
                    end
                    push!(queue, s)
                    push!(seen, s)
                end
            end
        end
    end
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

colors::Dict{Int,Union{Int,Symbol}} = Dict()

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
    do_color = haskey(colors, constraints(i))

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

translate_external::Bool = true
translate_internal::Bool = true

export set_color

"""
    set_color(s, color)

Sets the color of a index space `s` to the specified color.
Allowed colors are given by the
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

export disable_color

"""
    disable_color(s)

Disables the coloring for index space `s`.
"""
function disable_color(s::Int)
    delete!(colors, s)
    nothing
end

export enable_external_index_translation

"""
    enable_external_index_translation()

Enables external indices to be printed with names according to their
constraints.

!!! warning
    This can cause some confusion when different terms in the same
    expression contains different index constraints, as this will cause
    the same index to be printed with different names.

    ```jldoctest
    julia> using SpinAdaptedSecondQuantization

    julia> ex = E(1, 2) * occupied(1, 2) + E(1, 2) * virtual(1, 2)
    E_ab
    + E_ij
    julia> disable_external_index_translation()

    julia> ex # Now it is clearer that they are the same indices
    E_₁₂ C(₁∈v, ₂∈v)
    + E_₁₂ C(₁∈o, ₂∈o)
    julia> enable_external_index_translation()

    julia> ex # back to default behavour
    E_ab
    + E_ij
    ```

Enabled by default. Disabled using [`disable_external_index_translation`](@ref).
"""
function enable_external_index_translation()
    global translate_external = true
    nothing
end

export enable_internal_index_translation

"""
    enable_internal_index_translation()

Enables internal indices to be printed with names according to their
constraints.

# Example

```jldoctest
julia> using SpinAdaptedSecondQuantization

julia> ex = ∑(real_tensor("h", 1, 2) * E(1, 2) * electron(1, 2), 1:2)
∑_pq(h_pq E_pq)
julia> disable_internal_index_translation()

julia> ex
∑_₁₂(h_₁₂ E_₁₂) C(₁∈g, ₂∈g)
julia> ex *= E(1, 2) * occupied(2) * virtual(1) # Indices shifted up
∑_₃₄(h_₃₄ E_₃₄ E_ai) C(₃∈g, ₄∈g)
julia> enable_internal_index_translation()

julia> ex # Shift in indices is invisible when translated
∑_pq(h_pq E_pq E_ai)
```

Enabled by default.
"""
function enable_internal_index_translation()
    global translate_internal = true
    nothing
end

export disable_external_index_translation

"""
    disable_external_index_translation()

Disables external indices to be printed with names according to their
constraints. This can be useful as a debugging tool when dealing with
terms with different constraints on external indices in the same expression.

See [`enable_external_index_translation`](@ref) for example.

Enabled by default. Re-enable using [`enable_external_index_translation`](@ref)
"""
function disable_external_index_translation()
    global translate_external = false
    nothing
end

export disable_internal_index_translation

"""
    disable_internal_index_translation()

Disables internal indices to be printed with names according to their
constraints. This can be useful as a debugging tool if summation indices are
getting reordered/changed in weird ways.

See [`enable_internal_index_translation`](@ref) for example.

Enabled by default. Re-enable using [`enable_internal_index_translation`](@ref)
"""
function disable_internal_index_translation()
    global translate_internal = false
    nothing
end
