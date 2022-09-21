
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
