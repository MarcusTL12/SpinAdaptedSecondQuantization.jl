export e

"""
    SingletExcitationOperator

Type for representing e_pqrs as a single operator
"""
struct SingletDoubleExcitationOperator <: Operator
    p::Int
    q::Int
    r::Int
    s::Int

    # Constructor to respect symmetry
    function SingletDoubleExcitationOperator(p, q, r, s)
        if (r, s) < (p, q)
            new(r, s, p, q)
        else
            new(p, q, r, s)
        end
    end
end

function Base.show(io::IO,
    (
        e, constraints, translation
    )::Tuple{SingletDoubleExcitationOperator,Constraints,IndexTranslation})
    print(io, "e_")
    print_mo_index(io, constraints, translation, e.p, e.q, e.r, e.s)
end

function exchange_indices(e::SingletDoubleExcitationOperator, mapping)
    SingletDoubleExcitationOperator(
        exchange_index(e.p, mapping),
        exchange_index(e.q, mapping),
        exchange_index(e.r, mapping),
        exchange_index(e.s, mapping),
    )
end

function get_all_indices(e::SingletDoubleExcitationOperator)
    (e.p, e.q, e.r, e.s)
end

function Base.:(==)(a::SingletDoubleExcitationOperator, b::SingletDoubleExcitationOperator)
    (a.p, a.q, a.r, a.s) == (b.p, b.q, b.r, b.s)
end

function Base.isless(a::SingletDoubleExcitationOperator, b::SingletDoubleExcitationOperator)
    (a.p, a.q, a.r, a.s) < (b.p, b.q, b.r, b.s)
end

"""
    e(p, q, r, s)

    Constructs an expression containing a double excitation operator.
```
"""
e(p, q, r, s) = Expression(SingletDoubleExcitationOperator(p, q, r, s))

function convert_to_elementary_operators(o::SingletDoubleExcitationOperator)
    Expression(
        [
        (fermiondag(o.p, σ)*fermiondag(o.r, τ)*
         fermion(o.s, τ)*fermion(o.q, σ))[1]
        for σ in (α, β), τ in (α, β)
    ]
    )
end

function act_on_ket(op::SingletDoubleExcitationOperator)
    act_on_ket(eE(op.p, op.q, op.r, op.s))
end

function Base.adjoint(op::SingletDoubleExcitationOperator)
    SingletDoubleExcitationOperator(op.s, op.r, op.q, op.p)
end
