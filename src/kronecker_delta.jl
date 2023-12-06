export delta, δ

"""
    KroneckerDelta

Type representing a Kronecker delta of two or more MO-indices.

!!! note
    When constructing a `KroneckerDelta` it will return the integer 1 if
    the delta does not include two or more distinkt indices.
"""
struct KroneckerDelta
    indices::Vector{Int}

    function KroneckerDelta(indices)
        indices = unique!(sort(indices))

        if length(indices) <= 1
            1
        else
            new(indices)
        end
    end
end

function KroneckerDelta(indices...)
    KroneckerDelta(collect(indices))
end

function Base.show(io::IO,
    (
        d, constraints, translation
    )::Tuple{KroneckerDelta,Constraints,IndexTranslation})
    print(io, "δ_")

    for p in d.indices
        print_mo_index(io, constraints, translation, p)
    end
end

function print_latex(io::IO,
    (
        d, constraints, translation
    )::Tuple{KroneckerDelta,Constraints,IndexTranslation})
    print(io, "\\delta_{")

    for p in d.indices
        print_latex_mo_index(io, constraints, translation, p)
    end

    print(io, "}")
end

"""
    delta(indices...)

Construct an expression consisting of a single Kronecker delta of the given
indices.
"""
delta(indices...) = Expression(KroneckerDelta(indices...))

"""
    δ(indices...)

Unicode alias for [`delta`](@ref).
"""
δ(indices...) = delta(indices...)

function Base.isless(d1::KroneckerDelta, d2::KroneckerDelta)
    d1.indices < d2.indices
end

function Base.:(==)(d1::KroneckerDelta, d2::KroneckerDelta)
    d1.indices == d2.indices
end

"""
    exchange_indices(d::KroneckerDelta, mapping)

Returns a new `KroneckerDelta` with indices exchanged according to the `mapping`
"""
function exchange_indices(d::KroneckerDelta, mapping)
    KroneckerDelta([exchange_index(p, mapping) for p in d.indices])
end

"""
    compact_deltas(deltas::Vector{KroneckerDelta})

Returns a new reduced `Vector{KroneckerDelta}` where there exists only one
delta per group of equal indices.

Example:
```
δ_pq δ_qr δ_st -> δ_pqr δst
```
"""
function compact_deltas(deltas::Vector{KroneckerDelta})
    forest = DisjointSets{Int}()

    function add_index(i)
        if !haskey(forest.intmap, i)
            push!(forest, i)
        end
    end

    function add_index(r, i)
        if !haskey(forest.intmap, i)
            push!(forest, i)
        end

        union!(forest, r, i)
    end

    for d in deltas
        r, rest = Iterators.peel(d.indices)
        add_index(r)

        for i in rest
            add_index(r, i)
        end
    end

    group_dict = Dict{Int,Vector{Int}}()

    for i in collect(forest)
        r = find_root!(forest, i)
        if haskey(group_dict, r)
            push!(group_dict[r], i)
        else
            group_dict[r] = [i]
        end
    end

    new_deltas = KroneckerDelta[]

    for v in values(group_dict)
        d = KroneckerDelta(v)
        if d == 0
            return 0
        elseif d == 1
            @warn "Fusing deltas should not be able to produce 1!"
        else
            push!(new_deltas, d)
        end
    end

    sort!(new_deltas)
end
