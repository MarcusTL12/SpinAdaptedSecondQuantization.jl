export boson, bosondag

"""
    Boson Operators

The basic boson type operator.
"""
struct BosonOperator <: Operator
    dag::Bool
end

function Base.show(io::IO,
    (b, _, _)::Tuple{BosonOperator,Constraints,IndexTranslation})
    dag = b.dag ? '†' : '⁻'
    print(io, 'b', dag)
end

function print_latex(io::IO,
    (b, _, _)::Tuple{BosonOperator,Constraints,IndexTranslation})
    print(io, 'b')
    if b.dag
        print(io, "^\\dagger")
    end
end

function exchange_indices(b::BosonOperator, _)
    b
end

function get_all_indices(::BosonOperator)
    ()
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

function act_on_ket(op::BosonOperator)
    (op.dag) * Expression(op)
end
