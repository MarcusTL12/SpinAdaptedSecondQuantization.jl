# Implement commutation relations here

# To implement a commutation relation you implement the function
# `reductive_commutator(a, b)` between the two operator types.
# The function should return (∓1, [a, b]±), where the first integer should
# be the sign change of commutation (1 for commutator, -1 for anticommutator)

# When implementing a commutation relation between two different operator types
# only one order is required, as a generic function takes care of the other
# i.e. implement only [a, b]± and not also [b, a]±

# E_pq

function reductive_commutator(
    a::SingletExcitationOperator,
    b::SingletExcitationOperator
)
    p = a.p
    q = a.q
    r = b.p
    s = b.q

    (1, δ(q, r) * E(p, s) - δ(p, s) * E(r, q))
end

# a_p, a†_p

function reductive_commutator(a::FermionOperator, b::FermionOperator)
    (-1, δ(a.p, b.p) * (a.spin == b.spin) * (a.dag != b.dag))
end

function reductive_commutator(e::SingletExcitationOperator, a::FermionOperator)
    p = e.p
    q = e.q
    r = a.p

    (1, if a.dag
        δ(q, r) * fermiondag(p, a.spin)
    else
        -δ(p, r) * fermion(q, a.spin)
    end)
end

# e_pqrs

function reductive_commutator(
    a::SingletExcitationOperator,
    b::SingletDoubleExcitationOperator
)
    m, n = a.p, a.q
    p, q, r, s = b.p, b.q, b.r, b.s

    (1,
        δ(p, n) * e(m, q, r, s) - δ(m, q) * e(p, n, r, s) +
        δ(r, n) * e(p, q, m, s) - δ(s, q) * e(p, q, r, n)
    )
end

function reductive_commutator(
    a::Operator,
    b::SingletDoubleExcitationOperator
)
    t, u, v, w = b.p, b.q, b.r, b.s

    1, commutator(Expression(a), eE(t, u, v, w))
end

# T_pq

function reductive_commutator(t::TripletExcitationOperator,
    e::SingletExcitationOperator)

    p = t.p
    q = t.q

    r = e.p
    s = e.q

    (1, δ(r, q) * τ(p, s) - δ(p, s) * τ(r, q))
end

function reductive_commutator(t1::TripletExcitationOperator,
    t2::TripletExcitationOperator)

    p = t1.p
    q = t1.q

    r = t2.p
    s = t2.q

    (1, δ(r, q) * E(p, s) - δ(p, s) * E(r, q))
end

function reductive_commutator(e::TripletExcitationOperator, a::FermionOperator)
    p = e.p
    q = e.q
    r = a.p

    (1, if a.dag
        δ(q, r) * fermiondag(p, a.spin)
    else
        -δ(p, r) * fermion(q, a.spin)
    end * (a.spin == α ? 1 : -1))
end

# b, b†

function reductive_commutator(a::BosonOperator, b::BosonOperator)
    if a.dag && !b.dag
        return (1, Expression(-1))
    elseif !a.dag && b.dag
        return (1, Expression(1))
    else
        return (1, zero(Expression{Int64}))
    end
end

function reductive_commutator(::SingletExcitationOperator, ::BosonOperator)
    return (1, zero(Expression{Int64}))
end

function reductive_commutator(::SingletDoubleExcitationOperator, ::BosonOperator)
    return (1, zero(Expression{Int64}))
end

function reductive_commutator(::FermionOperator, ::BosonOperator)
    return (1, zero(Expression{Int64}))
end

function reductive_commutator(::TripletExcitationOperator, ::BosonOperator)
    return (1, zero(Expression{Int64}))
end
