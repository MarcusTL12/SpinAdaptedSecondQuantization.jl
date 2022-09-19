
"""
    SingletExcitationOperator

The basic E_pq type operator.
"""
struct SingletExcitationOperator <: Operator
    p::MOIndex
    q::MOIndex
end
