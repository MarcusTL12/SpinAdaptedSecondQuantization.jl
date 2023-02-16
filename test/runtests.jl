using Test

using SpinAdaptedSecondQuantization
import SpinAdaptedSecondQuantization as SASQ

@testset "indices" begin
    p = general(1)
    i = occupied(1)
    a = virtual(1)
    p1 = general(8)

    @test (space(p) <: GeneralOrbital) == true
    @test (space(p) <: OccupiedOrbital) == false
    @test (space(i) <: GeneralOrbital) == true
    @test (space(i) <: OccupiedOrbital) == true
 
    @test (space(i) == GeneralOrbital) == false
    @test (space(i) == OccupiedOrbital) == true
 
    @test isdisjoint(p, i) == false
    @test isdisjoint(a, i) == true
 
    @test isoccupied(i) == true
    @test isvirtual(i) == false
    @test isoccupied(a) == false
    @test isvirtual(a) == true
 
    @test string(p) == "p"
    @test string(p1) == "p₁"

    @test string(a) == "\e[36ma\e[39m"
    @test string(i) == "\e[92mi\e[39m"

    @test typeintersect(OccupiedOrbital, VirtualOrbital) == Union{}

end

@testset "kronecker delta" begin
    p = general(1)
    q = general(2)
    i = occupied(1)
    a = virtual(1)

    dpa = SASQ.KroneckerDelta(p, a)
    dia = SASQ.KroneckerDelta(i, a)
    dip = SASQ.KroneckerDelta(i, p)
    dpp = SASQ.KroneckerDelta(p, p)
    dpqi = SASQ.KroneckerDelta(p, q, i)

    @test string(dpa) == "δ_p\e[36ma\e[39m"
    @test string(dia) == "0"
    @test string(dip) == "δ_p\e[92mi\e[39m"
    @test string(dpp) == "1"
    @test string(dpqi) == "δ_pq\e[92mi\e[39m"

    @test dpa != 0
    @test dia == 0
    @test dip != 0
    @test dpp == 1
    @test dpqi != 0
end

@testset "kronecker delta compact" begin
    p = general(1)
    q = general(2)
    r = general(3)
    i = occupied(1)
    a = virtual(1)

    dpa = SASQ.KroneckerDelta(p, a)
    dri = SASQ.KroneckerDelta(r, i)
    dpq = SASQ.KroneckerDelta(p, q)
    dqr = SASQ.KroneckerDelta(q, r)

    v = [dpa, dpq, dri]

    @test string(v) == "SpinAdaptedSecondQuantization.KroneckerDelta[δ_p\e[36ma\e[39m, δ_pq, δ_r\e[92mi\e[39m]"

    v = SASQ.compact_deltas(v)

    @test string(v) == "SpinAdaptedSecondQuantization.KroneckerDelta[δ_pq\e[36ma\e[39m, δ_r\e[92mi\e[39m]"
    @test v != 0
    v = [dpa, dpq, dri, dqr]

    @test string(v) == "SpinAdaptedSecondQuantization.KroneckerDelta[δ_p\e[36ma\e[39m, δ_pq, δ_r\e[92mi\e[39m, δ_qr]"

    v = SASQ.compact_deltas(v)

    @test string(v) == "0"
    @test v == 0
end

@testset "singlet excitation operator" begin
    p = general(1)
    q = general(2)
    i = occupied(1)
    a = virtual(1)

    Epq = SASQ.SingletExcitationOperator(p, q)
    Eai = SASQ.SingletExcitationOperator(a, i)
    Eip = SASQ.SingletExcitationOperator(i, p)

    @test string(Epq) == "E_pq"
    @test string(Eai) == "E_\e[36ma\e[39m\e[92mi\e[39m"
    @test string(Eip) == "E_\e[92mi\e[39mp"
end

@testset "real tensor" begin
    p = general(1)
    q = general(2)

    hpq = SASQ.RealTensor("h", [p, q])

    @test string(hpq) == "h_pq"
end

@testset "term show" begin
    p = general(1)
    q = general(2)
    i = occupied(1)
    a = virtual(1)

    t = SASQ.Term(
        3 // 5,
        [i, a],
        [SASQ.KroneckerDelta(p, a),
            SASQ.KroneckerDelta(a, q)],
        [
            SASQ.RealTensor("h", [p, a]),
            SASQ.RealTensor("g", [i, a, i, i])
        ],
        [SASQ.SingletExcitationOperator(p, q)]
    )

    @test string(t) == "3/5 ∑_\e[92mi\e[39m\e[36ma\e[39m(δ_pq\e[36ma\e[39m g_\e[92mi\e[39m\e[36ma\e[39m\e[92mi\e[39m\e[92mi\e[39m h_p\e[36ma\e[39m E_pq)"
    @test string(SASQ.get_all_indices(t)) == "SpinAdaptedSecondQuantization.MOIndex[p, q, \e[92mi\e[39m, \e[36ma\e[39m]"

    t = SASQ.lower_delta_indices(t)

    @test string(t) == "3/5 ∑_\e[92mi\e[39m\e[36ma\e[39m(δ_pq\e[36ma\e[39m g_\e[92mi\e[39mp\e[92mi\e[39m\e[92mi\e[39m h_pp E_pp)"
end

@testset "expression constructors" begin
    p = general(1)
    q = general(2)
    i = occupied(1)
    a = virtual(1)

    Epq = E(p, q)
    dpi = δ(p, i)
    dai = δ(a, i)
    hai = real_tensor("h", a, i)

    @test string(Epq) == "E_pq"
    @test string(dpi) == "δ_p\e[92mi\e[39m"
    @test string(dai) == "0"
    @test string(hai) == "h_\e[36ma\e[39m\e[92mi\e[39m"
end

@testset "term exchange_indices" begin
    p = general(1)
    q = general(2)
    r = general(3)
    i = occupied(1)
    j = occupied(2)
    a = virtual(1)
    b = virtual(2)

    t = SASQ.Term(
        3 // 5,
        [i, a],
        [SASQ.KroneckerDelta(p, a)],
        [
            SASQ.RealTensor("h", [p, a]),
            SASQ.RealTensor("g", [i, a, i, q])
        ],
        [SASQ.SingletExcitationOperator(p, q)],
    )

    @test string(t) == "3/5 ∑_\e[92mi\e[39m\e[36ma\e[39m(δ_p\e[36ma\e[39m g_\e[92mi\e[39m\e[36ma\e[39m\e[92mi\e[39mq h_p\e[36ma\e[39m E_pq)"

    t2 = SASQ.exchange_indices(
        t,
        [p => q, q => p, a => b, i => j]
    )

    @test string(t2) == "3/5 ∑_\e[92mj\e[39m\e[36mb\e[39m(δ_q\e[36mb\e[39m \
g_\e[92mj\e[39m\e[36mb\e[39m\e[92mj\e[39mp h_q\e[36mb\e[39m E_qp)"

    t3 = SASQ.make_space_for_index(t, i)
    @test string(t3) == "3/5 ∑_\e[92mj\e[39m\e[36ma\e[39m(δ_p\e[36ma\e[39m g_\e[92mj\e[39m\e[36ma\e[39m\e[92mj\e[39mq h_p\e[36ma\e[39m E_pq)"

    t4 = SASQ.make_space_for_indices(t, [i, j, a])
    @test string(t4) == "3/5 ∑_\e[92mj\e[39m\e[36mb\e[39m(δ_p\e[36mb\e[39m g_\e[92mj\e[39m\e[36mb\e[39m\e[92mj\e[39mq h_p\e[36mb\e[39m E_pq)"

    t = SASQ.Term(
        3 // 7,
        SASQ.MOIndex[],
        [SASQ.KroneckerDelta(p, q, r)],
        [
            SASQ.RealTensor("h", [p, a]),
            SASQ.RealTensor("g", [i, a, i, q])
        ],
        [SASQ.SingletExcitationOperator(p, q)],
    )

    @test string(t) == "3/7 δ_pqr g_\e[92mi\e[39m\e[36ma\e[39m\e[92mi\e[39mq h_p\e[36ma\e[39m E_pq"

    t = SASQ.lower_delta_indices(t)

    @test string(t) == "3/7 δ_pqr g_\e[92mi\e[39m\e[36ma\e[39m\e[92mi\e[39mp h_p\e[36ma\e[39m E_pp"

    t2 = SASQ.exchange_indices(t, [p => q])

    @test string(t2) == "3/7 δ_qr g_\e[92mi\e[39m\e[36ma\e[39m\e[92mi\e[39mq h_q\e[36ma\e[39m E_qq"

    t3 = SASQ.exchange_indices(t2, [q => r])

    @test string(t3) == "3/7 g_\e[92mi\e[39m\e[36ma\e[39m\e[92mi\e[39mr h_r\e[36ma\e[39m E_rr"

    @test t4 == SASQ.exchange_indices(t4, [r => p])
end

@testset "term summation delta" begin
    p = general(1)
    q = general(2)
    i = occupied(1)
    j = occupied(2)
    a = virtual(1)
    b = virtual(2)

    t = SASQ.Term(
        1,
        SASQ.MOIndex[],
        [SASQ.KroneckerDelta(a, p)],
        [SASQ.RealTensor("h", [a, i])],
        SASQ.Operator[
            SASQ.SingletExcitationOperator(p, q)
        ],
    )

    t = SASQ.lower_delta_indices(t)

    #p -> a
    t1 = ∑(t, [p])
    @test t1 == ∑(δ(a,p) * real_tensor("h", p, i) * E(p, q), [p]).terms[1]
    t1 = SASQ.simplify_summation_deltas(t1)
    @test t1 == (real_tensor("h", a, i) * E(a, q)).terms[1]
    t2 = ∑(t1, [a])
    @test t2 == ∑(real_tensor("h", a, i) * E(a, q), [a]).terms[1]

    #a -> p
    t1 = ∑(t, [a])
    @test t1 == ∑(δ(a,p) * real_tensor("h", p, i) * E(p, q), [a]).terms[1]
    t1 = SASQ.simplify_summation_deltas(t1)
    t2 = simplify(∑(t1, [p]))
    @test t2 == ∑(real_tensor("h", a, i) * E(a, q), [a]).terms[1]

    #ap
    t1 = ∑(t, [a, p])
    @test t1 == ∑(δ(a,p) * real_tensor("h", p, i) * E(p, q), [p, a]).terms[1]
    t1 = SASQ.simplify_summation_deltas(t1)
    @test t1 == ∑(real_tensor("h", a, i) * E(a, q), [a]).terms[1]
end

@testset "term summation simplify" begin
    p = general(1)
    q = general(2)
    i = occupied(1)
    j = occupied(2)
    a = virtual(1)
    b = virtual(2)
    d = virtual(4)

    t = real_tensor("g", b, i, d, j)

    t = ∑(t, [b, d, i, j])
    t = SASQ.lower_summation_indices(t.terms[1])
    t = SASQ.sort_summation_indices(t)
    @test t == ∑(real_tensor("g", a, i, b, j), [i, j, a, b]).terms[1]
end

@testset "term exchange_indices constraints" begin
    p = general(1)
    q = general(2)
    r = general(3)
    s = general(4)
    i = occupied(1)
    j = occupied(2)
    a = virtual(1)
    b = virtual(2)

    t = SASQ.Term(
        3 // 5,
        SASQ.MOIndex[],
        SASQ.KroneckerDelta[],
        [
            SASQ.RealTensor("h", [a, b]),
        ],
        [SASQ.SingletExcitationOperator(i, j)],
    )

    t2 = SASQ.exchange_indices(
        t,
        [a => p, b => q, i => r, j => s]
    )

    @test t2 == (3//5 * real_tensor("h", p, q) * E(r, s) * constrain(p => VirtualOrbital, 
                                                                     q => VirtualOrbital,
                                                                     r => OccupiedOrbital,
                                                                     s => OccupiedOrbital)).terms[1]

    t3 = SASQ.exchange_indices(
        t2,
        [p => b]
    )

    @test t3 == (3//5 * real_tensor("h", b, q) * E(r, s) * constrain(q => VirtualOrbital,
                                                                     r => OccupiedOrbital,
                                                                     s => OccupiedOrbital)).terms[1]

    t4 = simplify(∑(t3, [s]))

    @test t4 == (∑(3//5 * real_tensor("h", b, q) * E(r, i) * constrain(q => VirtualOrbital,
                                                                       r => OccupiedOrbital), [i])).terms[1]
end

@testset "term multiplication" begin
    p = general(1)
    q = general(2)
    r = general(3)
    s = general(4)
    i = occupied(1)
    j = occupied(2)
    k = occupied(3)
    l = occupied(4)
    a = virtual(1)
    b = virtual(2)

    t1 = SASQ.Term(
        3 // 5,
        SASQ.MOIndex[i, j],
        SASQ.KroneckerDelta[],
        SASQ.Tensor[],
        [SASQ.SingletExcitationOperator(i, j)],
    )

    @test t1 * t1 == SASQ.Term(9//25, 
                         SASQ.MOIndex[i, j, k, l],
                         SASQ.KroneckerDelta[],
                         SASQ.Tensor[],
                         [SASQ.SingletExcitationOperator(i, j), SASQ.SingletExcitationOperator(k, l)],
                     )
end

@testset "expression addition" begin
    p = general(1)
    q = general(2)
    r = general(3)
    s = general(4)
    i = occupied(1)
    j = occupied(2)
    a = virtual(1)
    b = virtual(2)

    t1 = SASQ.Term(
        3 // 5,
        SASQ.MOIndex[],
        [
            SASQ.KroneckerDelta(i, p)
        ],
        SASQ.Tensor[],
        [SASQ.SingletExcitationOperator(p, q)],
    )

    e1 = δ(a, b)
    e2 = real_tensor("h", i, j)

    @test (e1 + 1) == (1 + e1)
    @test (e1 + e2) == (e2 + e1)
    @test (e1 + e2 + 2) == (2 + e1 + e2) == (e2 + 2 + e1)
end

@testset "expression multiplication" begin
    p = general(1)
    q = general(2)
    r = general(3)
    s = general(4)
    i = occupied(1)
    j = occupied(2)
    k = occupied(3)
    l = occupied(4)
    m = occupied(5)
    n = occupied(6)
    a = virtual(1)
    b = virtual(2)

    e1 = δ(a, b)
    e2 = real_tensor("h", i, j)
    h = ∑(real_tensor("h", i, j) * E(i, j), [i, j])
    g = ∑(real_tensor("g", i, j, k, l) * E(i, j) * E(k, l), [i, j, k, l])

    @test ((e1 + e2) * 3) // 5 * e2 == 3//5 * δ(a,b) * real_tensor("h", i, j) + 3//5 * real_tensor("h", i, j) * real_tensor("h", i, j)
    @test (h + h) == 2 * h
    @test simplify(g * h) == ∑(real_tensor("g", i, j, k, l) * real_tensor("h", m, n) * E(i, j) * E(k, l) * E(m, n), [i,j,k,l,m,n])
    @test simplify(h * g) == ∑(real_tensor("g", i, j, k, l) * real_tensor("h", m, n) * E(m, n) * E(i, j) * E(k, l), [i,j,k,l,m,n])
end

@testset "simple commutator" begin
    p = general(1)
    q = general(2)
    r = general(3)
    s = general(4)
    i = occupied(1)
    j = occupied(2)
    k = occupied(3)
    a = virtual(1)
    b = virtual(2)

    @test commutator(E(p, q), E(r, s)) == E(p, s) * δ(q, r) - E(r, q) * δ(p, s)
    @test commutator(E(i, j), E(a, k)) == - E(a, j) * δ(i, k)
    @test commutator(E(i, a), E(a, i)) == E(i, i) - E(a, a)
    @test commutator(E(a, i), E(a, i)) == SASQ.Expression(0)
    
    h = ∑(real_tensor("h", p, q) * E(p, q), [p, q])
    hcomm = simplify(commutator(h, E(p, q)))
    @test hcomm == summation(- real_tensor("h", q, r) * E(p, r) + real_tensor("h", r, p) * E(r, q), [r])
    @test hcomm * E(r, r) == summation(- real_tensor("h", q, s) * E(p, s) * E(r, r) + real_tensor("h", s, p) * E(s, q) * E(r, r), [s])
end

@testset "constrain" begin
    p = general(1)
    q = general(2)
    r = general(3)
    s = general(4)
    i = occupied(1)
    j = occupied(2)
    k = occupied(3)
    l = occupied(4)
    a = virtual(1)
    b = virtual(2)

    @test constrain(p => GeneralOrbital) == SASQ.Expression(1)
    @test string(constrain(p => OccupiedOrbital)) == "1 C(p∈O)"
    @test string(constrain(p => VirtualOrbital)) == "1 C(p∈V)"
    @test constrain(p => OccupiedOrbital) * constrain(p => VirtualOrbital) == SASQ.Expression(0)
    @test constrain(p => VirtualOrbital) * δ(p, a) == δ(p,a)
    @test constrain(p => VirtualOrbital) * δ(p, i) == SASQ.Expression(0)
    @test string(δ(p, q) * constrain(q => VirtualOrbital)) == "δ_pq C(p∈V)"
    @test (δ(p, q) * constrain(p => OccupiedOrbital)) * constrain(q => VirtualOrbital) == SASQ.Expression(0)
end

@testset "hf_expectation_value" begin
    p = general(1)
    q = general(2)
    r = general(3)
    s = general(4)
    i = occupied(1)
    j = occupied(2)
    k = occupied(3)
    l = occupied(4)
    a = virtual(1)
    b = virtual(2)

    @test hf_expectation_value(E(p, q)) == 2 * δ(p,q) * constrain(p => OccupiedOrbital)
    @test hf_expectation_value(E(i, j)) == 2 * δ(i,j)
    @test hf_expectation_value(E(a, b)) == SASQ.Expression(0)

    @test hf_expectation_value(E(i, j) * E(k, l)) == 4 * δ(i,j) * δ(k,l)
    @test hf_expectation_value(E(i, a) * E(b, j)) == 2 * δ(i,j) * δ(a,b)
    @test hf_expectation_value(E(p, q) * E(r, s)) == 4 * δ(p,q) * δ(r,s) * constrain(p=>OccupiedOrbital, r=>OccupiedOrbital) +
                                                     2 * δ(p,s) * δ(q,r) * constrain(p=>OccupiedOrbital, q=>VirtualOrbital)

    # Issue #13 from ExcitationOperators.jl
    e1 = hf_expectation_value(∑(E(k,j)*E(j,i), [i, j, k]))
    e2 = hf_expectation_value(∑(E(i,j)*E(j,i), [i, j])) 
    @test simplify(e1) == simplify(e2)
    
    @test hf_expectation_value(E(i,j)*E(j,i)) == 4 * δ(i,j)
end

# TODO add commented test-functions when simplify is improved and psym-integrals are added 
@testset "hf_equations" begin
    p = general(1)
    q = general(2)
    r = general(3)
    s = general(4)

    h = ∑(real_tensor("h", p, q) * E(p, q), [p, q])
    g = 1//2 * ∑(real_tensor("g", p, q, r, s) * e(p, q, r, s), [p, q, r, s])
    H = simplify(h + g)

    energy = simplify(hf_expectation_value(H))
    @show energy
    #@test energy == 2//1 * ∑(real_tensor(h, i, i), [i]) + ∑(2//1 * real_tensor("g", i, i, j, j) - real_tensor("g", i, j, j, i), [i, j])

    gradient = simplify(hf_expectation_value(commutator(H, E(p,q) - E(q,p))))
    #fock = real_tensor("h", p, q) + ∑(2//1 * real_tensor("g", p, q, i, i) - real_tensor("g", p, i, i, q), [i])

    @test gradient * constrain(p => OccupiedOrbital, q => OccupiedOrbital) == SASQ.Expression(0)
    @test gradient * constrain(p => VirtualOrbital, q => VirtualOrbital) == SASQ.Expression(0)
    #@test gradient * constrain(p => VirtualOrbital, q => OccupiedOrbital) == 4//1 * fock * constrain(p => VirtualOrbital, q => OccupiedOrbital)
end

@testset "jacobi_identity" begin
    p = general(1)
    q = general(2)
    r = general(3)
    s = general(4)
    t = general(5)
    u = general(6)
    
    f = commutator
    A = E(p,q)
    B = E(r,s)
    C = E(t,u)

    @test f(A, f(B, C)) + f(C, f(A, B)) + f(B, f(C, A)) == SASQ.Expression(0)
end
