export boson, bosondag

"""
    Boson Operators

The basic boson type operator.
"""

struct BosonOperator <: Operator
    dag::Bool
end

function Base.print(io::IO, constraints::Constraints, b::BosonOperator)
    dag = b.dag ? '†' : '⁻'
    print(io, 'b', dag)
end

function exchange_indices(b::BosonOperator, mapping)
    b
end

function get_all_indices(b::BosonOperator)
    []
end

function Base.isless(a::BosonOperator, b::BosonOperator)
    a.dag < b.dag
end

function Base.:(==)(a::BosonOperator, b::BosonOperator)
    a.dag == b.dag
end

# Externally visible constructor
boson() = Expression(BosonOperator(false))
bosondag() = Expression(BosonOperator(true))

function Base.:(==)(a::BosonOperator, b::SingletExcitationOperator)
    false
end

function act_on_ket(op::BosonOperator)
    (op.dag) * Expression(op)
end
