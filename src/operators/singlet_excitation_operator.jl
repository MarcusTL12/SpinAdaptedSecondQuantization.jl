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

# Externally visible constructor
E(p, q) = Term(
    1,
    MOIndex[],
    KroneckerDelta[],
    Tensor[],
    Operator[SingletExcitationOperator(p, q)]
)
