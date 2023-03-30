export fermion, fermiondag

"""
    Fermionic Operators

The basic fermionic type operator.
"""

struct FermionOperator <: Operator
    p::Int
    spin::Spin
    dag::Bool
end

function Base.print(io::IO, constraints::Constraints,
    translation::IndexTranslation, a::FermionOperator)
    dag = a.dag ? '†' : '⁻'
    print(io, 'a', dag, '_')
    print_mo_index(io, constraints, translation, a.p)
    print(io, a.spin)
end

function exchange_indices(a::FermionOperator, mapping)
    FermionOperator(exchange_index(a.p, mapping), a.spin, a.dag)
end

function get_all_indices(e::FermionOperator)
    (e.p,)
end

function Base.isless(a::FermionOperator, b::FermionOperator)
    (b.dag, a.spin, a.p) < (a.dag, b.spin, b.p)
end

function Base.:(==)(a::FermionOperator, b::FermionOperator)
    (a.p, a.spin, a.dag) == (b.p, b.spin, b.dag)
end

# Externally visible constructor
fermion(p, spin) = Expression(FermionOperator(p, spin, false))
fermiondag(p, spin) = Expression(FermionOperator(p, spin, true))

function act_on_ket(op::FermionOperator)
    Expression(op) * if op.dag
        virtual(op.p)
    else
        occupied(op.p)
    end
end
