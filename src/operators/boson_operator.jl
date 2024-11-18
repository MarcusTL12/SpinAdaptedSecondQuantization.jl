export boson, bosondag

"""
    Boson Operators

The basic boson type operator.
"""
struct BosonOperator <: Operator
    p::Int
    dag::Bool
end

function Base.show(io::IO,
    (
        b, constraints, translation
    )::Tuple{BosonOperator,Constraints,IndexTranslation})
    dag = b.dag ? '†' : '⁻'
    print(io, 'b', dag, '_')
    print_mo_index(io, constraints, translation, b.p)
end

function exchange_indices(b::BosonOperator, mapping)
    BosonOperator(exchange_index(b.p, mapping), b.dag)
end

function get_all_indices(b::BosonOperator)
    (b.p,)
end

function Base.isless(a::BosonOperator, b::BosonOperator)
    (a.dag, a.p) < (b.dag, b.p)
end

function Base.:(==)(a::BosonOperator, b::BosonOperator)
    (a.dag, a.p) == (b.dag, b.p)
end

# Externally visible constructor
boson(p) = Expression(BosonOperator(p, false))
bosondag(p) = Expression(BosonOperator(p, true))

function act_on_ket(op::BosonOperator)
    (op.dag) * Expression(op)
end

function Base.adjoint(op::BosonOperator)
    BosonOperator(op.p, !op.dag)
end
