export fermion, fermiondag

"""
    Fermionic Operators

The basic fermionic type operator.
"""

struct FermionOperator <: Operator
    p::Int
    spin::Bool
    dag::Bool
end

function Base.print(io::IO, constraints::Constraints, a::FermionOperator)
    spin = a.spin ? 'α' : 'β'
    dag = a.dag ? '†' : '⁻'
    print(io, 'a', dag, '_')
    print_mo_index(io, constraints, a.p)
    print(io, spin)
end

function exchange_indices(a::FermionOperator, mapping)
    FermionOperator(exchange_index(a.p, mapping), a.spin, a.dag)
end

function get_all_indices(e::FermionOperator)
    [e.p]
end

function Base.isless(a::FermionOperator, b::FermionOperator)
    (a.dag, a.spin, a.p) < (b.dag, b.spin, b.p)
end

function Base.:(==)(a::FermionOperator, b::FermionOperator)
    (a.p, a.spin, a.dag) == (b.p, b.spin, b.dag)
end

# Externally visible constructor
fermion(p, spin) = Expression(FermionOperator(p, spin, false))
fermiondag(p, spin) = Expression(FermionOperator(p, spin, true))

# Commutation relations:
function commutator(a::FermionOperator, b::FermionOperator)
    anticommutator(a, b) - 2 * Expression(b) * Expression(a)
end

function anticommutator(a::FermionOperator, b::FermionOperator)
    δ(a.p, b.p) * (a.spin == b.spin) * (a.dag != b.dag)
end
