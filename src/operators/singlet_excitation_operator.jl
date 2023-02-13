export E, e

"""
    SingletExcitationOperator

The basic E_pq type operator.
"""
struct SingletExcitationOperator <: Operator
    p::MOIndex
    q::MOIndex
end

function Base.show(io::IO, e::SingletExcitationOperator)
    print(io, "E_", e.p, e.q)
end

function exchange_indices(e::SingletExcitationOperator, mapping)
    SingletExcitationOperator(
        exchange_index(e.p, mapping),
        exchange_index(e.q, mapping)
    )
end

function get_all_indices(e::SingletExcitationOperator)
    [e.p, e.q]
end

# Externally visible constructor
E(p, q) = Expression(SingletExcitationOperator(p, q))

e(p, q, r, s) = E(p, q) * E(r, s) - δ(r, q) * E(p, s)

# Commutation relations:

function commutator(a::SingletExcitationOperator, b::SingletExcitationOperator)
    p = a.p
    q = a.q
    r = b.p
    s = b.q

    δ(q, r) * E(p, s) - δ(p, s) * E(r, q)
end
