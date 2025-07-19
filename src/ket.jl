export act_on_ket

"""
    act_on_ket(ex::Expression, [max_ops=Inf])

Project `ex` onto a ``|\\text{HF}\\rangle`` ket to the right.
This is done termwise and for each term the folliwing recursive formula is used
to project an arbitrary string of operators on a ket.

The base case is that there is an implementation of `act_on_ket(o)` for each
operator type in the string such that we have

``
A_i |\\text{HF}\\rangle = \\tilde A_i |\\text{HF}\\rangle
``

where ``\\tilde A_i`` is constrained and simplified from the projection.
The recursive step is then

``
A_1 A_2 ... A_n |\\text{HF}\\rangle
=
A_1 A_2 ... \\tilde A_n |\\text{HF}\\rangle
=
[A_1 A_2 ..., \\tilde A_n]_\\Gamma |\\text{HF}\\rangle +
\\Gamma \\tilde A_n (A_1 A_2 ... A_{n-1}) |\\text{HF}\\rangle
``

where the ``[A, B]_\\Gamma`` notation is the [`reductive_commutator`](@ref)
which chooses whether to compute a commutatior or anticommutator depending on
the input to reduce the number of operator in the output.

The optional parameter `max_ops` specifies the maximum number of operator to
keep in any terms. If only a few excitations are needed, specifying this can
greatly speed up the projection.

# Example

```jldoctest
julia> using SpinAdaptedSecondQuantization

julia> E(1, 2) * E(3, 4) * electron(1:4...)
E_pq E_rs
julia> ket = simplify(act_on_ket(ans))
E_ai E_bj
+ 2 δ_ij E_ak
- δ_ik E_aj
+ δ_bc E_ai
+ 2 δ_jk E_ai
+ 4 δ_ij δ_kl
+ 2 δ_ij δ_ab
julia> disable_external_index_translation()

julia> ket
E_₁₂ E_₃₄ C(₁∈v, ₂∈o, ₃∈v, ₄∈o)
+ 2 δ_₁₂ E_₃₄ C(₁∈o, ₂∈o, ₃∈v, ₄∈o)
- δ_₁₄ E_₃₂ C(₁∈o, ₂∈o, ₃∈v, ₄∈o)
+ δ_₂₃ E_₁₄ C(₁∈v, ₂∈v, ₃∈v, ₄∈o)
+ 2 δ_₃₄ E_₁₂ C(₁∈v, ₂∈o, ₃∈o, ₄∈o)
+ 4 δ_₁₂ δ_₃₄ C(₁∈o, ₂∈o, ₃∈o, ₄∈o)
+ 2 δ_₁₄ δ_₂₃ C(₁∈o, ₂∈v, ₃∈v, ₄∈o)
julia> enable_external_index_translation()
```
"""
function act_on_ket(ex::Expression{T}, max_ops=Inf) where {T}
    nth = Threads.nthreads()
    terms = [Term{T}[] for _ in 1:nth]
    Threads.@threads for id in 1:nth
        for i in id:nth:length(ex.terms)
            append!(terms[id], act_on_ket(ex[i], max_ops).terms)
        end
    end

    all_terms, rest = Iterators.peel(terms)
    for other_terms in rest
        append!(all_terms, other_terms)
    end

    Expression(all_terms)
end

function act_on_ket_unthreaded(ex::Expression{T}, max_ops=Inf) where {T}
    nth = Threads.nthreads()
    terms = [Term{T}[] for _ in 1:nth]
    for id in 1:nth
        for i in id:nth:length(ex.terms)
            append!(terms[id], act_on_ket(ex[i], max_ops).terms)
        end
    end

    all_terms, rest = Iterators.peel(terms)
    for other_terms in rest
        append!(all_terms, other_terms)
    end

    Expression(all_terms)
end

function act_on_ket(t::Term{A}, max_ops) where {A<:Number}
    if iszero(t.scalar)
        return Expression(zero(A))
    end
    if isempty(t.operators)
        return Expression([t])
    end

    copyt = copy(t)
    right_op = pop!(copyt.operators)
    right_op_act = act_on_ket(right_op)
    copyt_act = act_on_ket(copyt,
        max_ops - minimum(length(t.operators) for t in right_op_act.terms))

    terms = Term{A}[]
    for r in right_op_act.terms
        Γ, comm = reductive_commutator_fuse(copyt, r)

        if length(r.operators) <= max_ops
            new_max = max_ops - length(r.operators)
            append!(terms, Γ * fuse(r, ter)
                           for ter in copyt_act.terms
                           if length(ter.operators) <= new_max)
        end

        append!(terms, act_on_ket_unthreaded(comm, max_ops).terms)
    end

    Expression(terms)
end

export hf_expectation_value

"""
    hf_expectation_value(ex::Expression) = act_on_ket(ex, 0)

Alias for calling [`act_on_ket`](@ref) keeping no operators, which gives the
HF expectation value.
"""
hf_expectation_value(ex::Expression) = act_on_ket(ex, 0)
