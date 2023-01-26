export E

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

# Externally visible constructor
E(p, q) = Expression(SingletExcitationOperator(p, q))
