export Q

struct PairFermionOperator <: Operator
    p::Int
    q::Int
    dag::Bool
    function PairFermionOperator(p, q, dag)
        new(minmax(p, q)..., dag)
    end
end

function Base.show(io::IO,
    (
        o, constraints, translation
    )::Tuple{PairFermionOperator,Constraints,IndexTranslation})
    dag = o.dag ? '†' : '⁻'
    print(io, 'Q', dag, '_')
    print_mo_index(io, constraints, translation, o.p)
    print_mo_index(io, constraints, translation, o.q)
end

function exchange_indices(o::PairFermionOperator, mapping)
    PairFermionOperator(exchange_index(o.p, mapping), exchange_index(o.q, mapping), o.dag)
end

function get_all_indices(o::PairFermionOperator)
    (o.p, o.q)
end

function Base.isless(a::PairFermionOperator, b::PairFermionOperator)
    (b.dag, a.p, a.q) < (a.dag, b.p, b.q)
end

function Base.:(==)(a::PairFermionOperator, b::PairFermionOperator)
    (a.p, a.q, a.dag) == (b.p, b.q, b.dag)
end

function Base.adjoint(op::PairFermionOperator)
    PairFermionOperator(op.p, op.q, !op.dag)
end

Q(p, q) = Expression(PairFermionOperator(p, q, false))

function act_on_ket(op::PairFermionOperator)
    Expression(op) * if op.dag
        virtual(op.p, op.q)
    else
        occupied(op.p, op.q)
    end
end
