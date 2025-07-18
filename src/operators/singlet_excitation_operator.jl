export E, e

"""
    SingletExcitationOperator

The basic E_pq type operator.
"""
struct SingletExcitationOperator <: Operator
    p::Int
    q::Int
end

function Base.show(io::IO,
    (
        e, constraints, translation
    )::Tuple{SingletExcitationOperator,Constraints,IndexTranslation})
    print(io, "E_")
    print_mo_index(io, constraints, translation, e.p, e.q)
end

function exchange_indices(e::SingletExcitationOperator, mapping)
    SingletExcitationOperator(
        exchange_index(e.p, mapping),
        exchange_index(e.q, mapping)
    )
end

function get_all_indices(e::SingletExcitationOperator)
    (e.p, e.q)
end

function Base.:(==)(a::SingletExcitationOperator, b::SingletExcitationOperator)
    (a.p, a.q) == (b.p, b.q)
end

function Base.isless(a::SingletExcitationOperator, b::SingletExcitationOperator)
    (a.p, a.q) < (b.p, b.q)
end

"""
    E(p, q)

Constructs an expression containing a single excitation operator.

# Example

```jldoctest
julia> using SpinAdaptedSecondQuantization

julia> E(1, 2) * electron(1, 2)
E_pq
julia> E(1, 2) * E(3, 4) * occupied(2, 4) * virtual(1, 3)
E_ai E_bj
```
"""
E(p, q) = Expression(SingletExcitationOperator(p, q))

"""
    e(p, q, r, s) = E(p, q) * E(r, s) - δ(r, q) * E(p, s)

Alias for the two electron singlet excitation operator.

# Example:

```jldoctest
julia> using SpinAdaptedSecondQuantization

julia> e(1, 2, 3, 4) * electron(1, 2, 3, 4)
E_pq E_rs
- δ_qr E_ps
```
"""
e(p, q, r, s) = E(p, q) * E(r, s) - δ(r, q) * E(p, s)

function convert_to_elementary_operators(o::SingletExcitationOperator)
    Expression(
        [(fermiondag(o.p, spin)*fermion(o.q, spin))[1] for spin in (α, β)]
    )
end

function act_on_ket(op::SingletExcitationOperator)
    p = op.p
    q = op.q
    E(p, q) * virtual(p) * occupied(q) +
    2 * δ(p, q) * occupied(p, q)
end

function Base.adjoint(op::SingletExcitationOperator)
    SingletExcitationOperator(op.q, op.p)
end
