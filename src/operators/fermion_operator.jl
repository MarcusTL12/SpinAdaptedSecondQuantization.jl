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

function Base.show(io::IO,
    (
        a, constraints, translation
    )::Tuple{FermionOperator,Constraints,IndexTranslation})
    dag = a.dag ? '†' : '⁻'
    print(io, 'a', dag, '_')
    print_mo_index(io, constraints, translation, a.p)
    print(io, a.spin)
end

function print_latex(io::IO,
    (
        a, constraints, translation
    )::Tuple{FermionOperator,Constraints,IndexTranslation})
    print(io, 'a')
    if a.dag
        print(io, "^\\dagger")
    end
    print(io, "_{")
    print_latex_mo_index(io, constraints, translation, a.p)
    if a.spin == α
        print(io, "\\alpha}")
    elseif a.spin == β
        print(io, "\\beta}")
    end
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
