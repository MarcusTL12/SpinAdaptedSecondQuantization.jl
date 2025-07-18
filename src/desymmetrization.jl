export symmetrize, desymmetrize, make_permutation_mappings

"""
    make_permutation_mappings(index_groups)

Convenience function to make mapping for all permutations of groups of indices.
See [`symmetrize`](@ref) and [`desymmetrize`](@ref).
"""
function make_permutation_mappings(index_groups)
    mappings = Vector{Pair{Int,Int}}[]

    for perm in PermGen(length(index_groups))
        pd = perm.data
        mapping = Pair{Int,Int}[]

        for (i, j) in enumerate(pd)
            for (from, to) in zip(index_groups[i], index_groups[j])
                push!(mapping, from => to)
            end
        end

        push!(mappings, mapping)
    end

    mappings
end

"""
    desymmetrize(ex::Expression, mappings)

Function to look for terms which are equal under the index permutations
given by the `mappings`. See [`symmetrize`](@ref) and
[`make_permutation_mappings`](@ref).

Returns three expressions:
1. The first expression contains the terms where redundant copies were found
    which are equal under the given symmetry. These needs to be symmetrized
    to obtain the original expression, usally after numerically evaluating
    the non-redundant terms.
2. The second expression constins the "self-symmetric" terms, which on their own
    are invariant under the given symmetries and do not need to be symmetrized.
3. The last expression contains any left over "non-symmetric" terms. This is
    expected to be zero when desymmetrizing an expression with a known symmetry
    such as CCSD omega equations.

# Examples

Setting up an example expression:
```jldoctest desymmetrize
julia> using SpinAdaptedSecondQuantization

julia> ∑(real_tensor("F", 1, 5) * psym_tensor("t", 3, 4, 5, 2) * \
occupied(2, 4) * virtual(1, 3, 5), [5])
∑_c(F_ac t_bjci)
julia> symmetrize(ans, make_permutation_mappings([(1, 2), (3, 4)]))
∑_c(F_ac t_bjci)
+ ∑_c(F_bc t_aicj)
julia> ans + ∑(psym_tensor("g", 1, 5, 3, 6) * psym_tensor("t", 5, 2, 6, 4) * \
occupied(2, 4) * virtual(1, 3, 5, 6), 5:6)
∑_c(F_ac t_bjci)
+ ∑_c(F_bc t_aicj)
+ ∑_cd(g_acbd t_cidj)
julia> ans + real_tensor("x", 1, 2, 3, 4) * occupied(2, 4) * virtual(1, 3)
x_aibj
+ ∑_c(F_ac t_bjci)
+ ∑_c(F_bc t_aicj)
+ ∑_cd(g_acbd t_cidj)
```

# Desymmetrizing:
# ```jldoctest desymmetrize
# julia> r, ss, ns = desymmetrize(ans * 1//1, \
# make_permutation_mappings([(1, 2), (3, 4)]));

# julia> r # Only one of the two redundant terms survives
# ∑_c(F_ac t_bjci)
# julia> ss # The self-symmetric term
# ∑_cd(g_acbd t_cidj)
# julia> ns # The non-symmetric term
# x_aibj
# ```

!!! warning
    `desymmetrize` does repeated calls to [`simplify_heavy`](@ref) which can
    be slow if terms have many (>=8) summation indices.
"""
function desymmetrize(ex_::Expression{T}, mappings) where {T<:Number}
    accounted_for = Set{Int}()

    new_T = promote_type(T, Rational{Int})

    ex = promote_scalar(new_T, ex_)

    non_symmetric = Term{new_T}[]
    self_symmetric = Term{new_T}[]
    symmetrize = Term{new_T}[]

    for (i, t) in enumerate(ex.terms)
        if i ∈ accounted_for
            continue
        end

        is_self_symmetric = true
        has_redundant_terms = false

        other_inds = Dict{Int,Int}()
        for mapping in mappings
            other_term = simplify_heavy(exchange_indices(t, mapping))
            if other_term == t
                other_inds[i] = get(other_inds, i, 0) + 1
            else
                is_self_symmetric = false
                for (j, t2) in enumerate(ex.terms)
                    if other_term == t2
                        has_redundant_terms = true
                        other_inds[j] = get(other_inds, j, 0) + 1
                        break
                    end
                end
            end
        end

        if is_self_symmetric
            push!(self_symmetric, t)
            push!(accounted_for, i)
        elseif has_redundant_terms
            @assert allequal(values(other_inds))
            push!(symmetrize, t * (1 // first(values(other_inds))))
            push!(accounted_for, i)
            union!(accounted_for, keys(other_inds))
        else
            @show Expression([t]) other_inds
            push!(non_symmetric, t)
            push!(accounted_for, i)
        end
    end

    non_symmetric = Expression(non_symmetric)
    self_symmetric = Expression(self_symmetric)
    symmetrize = Expression(symmetrize)

    symmetrize, self_symmetric, non_symmetric
end

"""
    symmetrize(ex::Expression, mappings)

Function to expand the expression `ex` including all permutations of indices
given by `mappings`. This is common to use to expand all permuations among
neighbouring pairs of indices, such as to compute
``P_{ijk}^{abc} \\Omega_{ijk}^{abc}``. For this specific purpose we provide the
convenience function [`make_permutation_mappings`](@ref) which produces all
permutations between groups of indices.

```@meta
DocTestSetup = :(using SpinAdaptedSecondQuantization)
```

```jldoctest
julia> x = real_tensor("t", 1, 2, 3, 4, 5, 6) *
           occupied(2, 4, 6) * virtual(1, 3, 5)
t_aibjck
julia> symmetrize(x, make_permutation_mappings([(1, 2), (3, 4), (5, 6)]))
t_aibjck
+ t_aickbj
+ t_bjaick
+ t_bjckai
+ t_ckaibj
+ t_ckbjai
```
"""
function symmetrize(ex::Expression{T}, mappings) where {T<:Number}
    terms = Term{T}[]

    for t in ex.terms
        for mapping in mappings
            other_term = exchange_indices(t, mapping)
            push!(terms, other_term)
        end
    end

    Expression(terms)
end
