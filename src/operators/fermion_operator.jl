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

function Base.show(io::IO, a::FermionOperator)
    spin = a.spin ? "↑" : "↓"
    dag = a.dag ? "+" : "-"
    print(io, "($(spin)$(dag))_")
    print_mo_index(io, a.p)
end

function exchange_indices(a::FermionOperator, mapping)
    FermionOperator(exchange_index(a.p, mapping), a.spin, a.dag)
end

function get_all_indices(e::FermionOperator)
    [e.p]
end

# Externally visible constructor
fermion(p, spin) = Expression(FermionOperator(p, spin, false))
fermiondag(p, spin) = Expression(FermionOperator(p, spin, true))

# Commutation relations:
function commutator(a::FermionOperator, b::FermionOperator)
    anticommutator(a, b) - 2 * Expression(b) * Expression(a)
end

function anticommutator(a::FermionOperator, b::FermionOperator)
    δ(a.p, b.p) * (a.spin==b.spin) * (a.dag != b.dag)
end