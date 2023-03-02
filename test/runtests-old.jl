using Test

using SpinAdaptedSecondQuantization
import SpinAdaptedSecondQuantization as SASQ

@testset "ortital spaces" begin
    @test OccupiedOrbital <: GeneralOrbital
    @test VirtualOrbital <: GeneralOrbital
    @test isdisjoint(VirtualOrbital, OccupiedOrbital)
    @test !isdisjoint(GeneralOrbital, OccupiedOrbital)
    @test typeintersect(GeneralOrbital, OccupiedOrbital) == OccupiedOrbital
    @test typeintersect(OccupiedOrbital, VirtualOrbital) == Union{}
    @test SpinAdaptedSecondQuantization.is_strict_subspace(OccupiedOrbital, GeneralOrbital)
    @test !SpinAdaptedSecondQuantization.is_strict_subspace(GeneralOrbital, GeneralOrbital)
end

@testset "print indices" begin
    p = general(1)
    p1 = general(8)
    i = occupied(1)
    a = virtual(1)
    @test string(p) == "p"
    @test string(p1) == "p₁"
    @test string(i) == "i"
    @test string(a) == "a"
end

@testset "kronecker delta" begin
    p = general(1)
    q = general(2)
    dpp = SASQ.KroneckerDelta(p, p)
    dpq = SASQ.KroneckerDelta(p, q)

    @test string(dpp) == "1"
    @test string(dpq) == "δ_pq"
    @test dpp == 1
end

@testset "kronecker delta compact" begin
    p = general(1)
    q = general(2)
    r = general(3)
    dpq = SASQ.KroneckerDelta(p, q)
    dqr = SASQ.KroneckerDelta(q, r)

    v = [dpq]
    @test string(v) == "SpinAdaptedSecondQuantization.KroneckerDelta[δ_pq]"

    v = [dpq, dqr]
    @test string(v) == "SpinAdaptedSecondQuantization.KroneckerDelta[δ_pq, δ_qr]"

    v = SASQ.compact_deltas(v)
    @test string(v) == "SpinAdaptedSecondQuantization.KroneckerDelta[δ_pqr]"
end

@testset "print singlet excitation operator" begin
    p = general(1)
    q = general(2)

    Epq = SASQ.SingletExcitationOperator(p, q)
    @test string(Epq) == "E_pq"
end

@testset "print real tensor" begin
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
        [SASQ.SingletExcitationOperator(p, q)],
        SASQ.Constraints(i => OccupiedOrbital, a => VirtualOrbital)
    )

    @test string(t) == "3/5 ∑_ia(δ_pqa g_iaii h_pa E_pq) C(p∈V, i∈O)"
    @test string(SASQ.get_all_indices(t)) == "SpinAdaptedSecondQuantization.MOIndex[p, q, i, a]"

    t = SASQ.lower_delta_indices(t)
    @test string(t) == "3/5 ∑_ia(δ_pqa g_ipii h_pp E_pp) C(p∈V, i∈O)"
end

@testset "term exchange_indices" begin
    p = general(1)
    q = general(2)
    r = general(3)
    i = occupied(1)
    a = virtual(1)

    t = SASQ.Term(
        3 // 5,
        [i, a],
        [SASQ.KroneckerDelta(p, a)],
        [
            SASQ.RealTensor("h", [p, a]),
            SASQ.RealTensor("g", [i, a, i, q])
        ],
        [SASQ.SingletExcitationOperator(p, q)],
        SASQ.Constraints(i => OccupiedOrbital, a => VirtualOrbital)
    )

    @test string(t) == "3/5 ∑_ia(δ_pa g_iaiq h_pa E_pq) C(p∈V, i∈O)"

    t = SASQ.Term(
        3 // 7,
        SASQ.MOIndex[],
        [SASQ.KroneckerDelta(p, q, r)],
        [
            SASQ.RealTensor("h", [p, a]),
            SASQ.RealTensor("g", [i, a, i, q])
        ],
        [SASQ.SingletExcitationOperator(p, q)],
        SASQ.Constraints(i => OccupiedOrbital, a => VirtualOrbital)
    )

    @test string(t) == "3/7 δ_pqr g_iaiq h_pa E_pq C(i∈O, a∈V)"

    t = SASQ.lower_delta_indices(t)

    @test string(t) == "3/7 δ_pqr g_iaip h_pa E_pp C(i∈O, a∈V)"

    t2 = SASQ.exchange_indices(t, [p => q])

    @test string(t2) == "3/7 δ_qr g_iaiq h_qa E_qq C(i∈O, a∈V)"

    t3 = SASQ.exchange_indices(t2, [q => r])

    @test string(t3) == "3/7 g_iair h_ra E_rr C(i∈O, a∈V)"
end

@testset "term summation delta" begin
    p = general(1)
    q = general(2)
    i = occupied(1)
    a = virtual(1)

    t = SASQ.Term(
        1,
        SASQ.MOIndex[],
        [SASQ.KroneckerDelta(a, p)],
        [SASQ.RealTensor("h", [a, i])],
        SASQ.Operator[
            SASQ.SingletExcitationOperator(p, q)
        ],
        SASQ.Constraints(i => OccupiedOrbital, a => VirtualOrbital)
    )

    t = SASQ.lower_delta_indices(t)

    #p -> a
    t1 = ∑(t, [p])
    t1 = SASQ.simplify_summation_deltas(t1)
    t1 = ∑(t1, [a]) * constrain(a => VirtualOrbital)[1]

    #a -> p
    t2 = ∑(t, [a]) * constrain(a => VirtualOrbital)[1]
    t2 = SASQ.simplify_summation_deltas(t2)
    t2 = simplify(∑(t2, [p]))

    #ap
    t3 = ∑(t, [a, p]) * constrain(a => VirtualOrbital)[1]
    t3 = simplify(t3)

    @test t1 == t2 == t3
end

@testset "term multiplication" begin
    p = general(1)
    q = general(2)
    r = general(3)
    s = general(4)

    t1 = SASQ.Term(
        3 // 5,
        SASQ.MOIndex[p, q],
        SASQ.KroneckerDelta[],
        SASQ.Tensor[],
        [SASQ.SingletExcitationOperator(p, q)],
    )

    @test simplify(t1 * t1) == SASQ.Term(9 // 25,
        SASQ.MOIndex[p, q, r, s],
        SASQ.KroneckerDelta[],
        SASQ.Tensor[],
        [SASQ.SingletExcitationOperator(p, q), SASQ.SingletExcitationOperator(r, s)],
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
    i = general(1)
    j = general(2)
    k = general(3)
    l = general(4)
    m = general(5)
    n = general(6)
    a = virtual(1)
    b = virtual(2)

    e1 = δ(a, b)
    e2 = real_tensor("h", i, j)
    h = ∑(real_tensor("h", i, j) * E(i, j), [i, j])
    g = ∑(real_tensor("g", i, j, k, l) * E(i, j) * E(k, l), [i, j, k, l])

    @test ((e1 + e2) * 3) // 5 * e2 == 3 // 5 * δ(a, b) * real_tensor("h", i, j) + 3 // 5 * real_tensor("h", i, j) * real_tensor("h", i, j)
    @test (h + h) == 2 * h
    @test simplify(g * h) == ∑(real_tensor("g", i, j, k, l) * real_tensor("h", m, n) * E(i, j) * E(k, l) * E(m, n), [i, j, k, l, m, n])
    @test simplify(h * g) == ∑(real_tensor("g", i, j, k, l) * real_tensor("h", m, n) * E(m, n) * E(i, j) * E(k, l), [i, j, k, l, m, n])
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
    @test commutator(E(i, a), E(a, i)) == E(i, i) - E(a, a)
    @test commutator(E(a, i), E(a, i)) == SASQ.Expression(0)

    h = ∑(real_tensor("h", p, q) * E(p, q), [p, q])
    hcomm = simplify(commutator(h, E(p, q)))
    @test hcomm == summation(-real_tensor("h", q, r) * E(p, r) + real_tensor("h", r, p) * E(r, q), [r])
    @test hcomm * E(r, r) == summation(-real_tensor("h", q, s) * E(p, s) * E(r, r) + real_tensor("h", s, p) * E(s, q) * E(r, r), [s])
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
    @test constrain(p => VirtualOrbital, a => VirtualOrbital) * δ(p, a) == δ(p, a) * constrain(p => VirtualOrbital)
    @test constrain(p => VirtualOrbital, i => OccupiedOrbital) * δ(p, i) == SASQ.Expression(0)
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

    @test hf_expectation_value(E(p, q)) == 2 * δ(p, q) * constrain(p => OccupiedOrbital)
    @test hf_expectation_value(E(p, q) * E(r, s)) == 4 * δ(p, q) * δ(r, s) * constrain(p => OccupiedOrbital, r => OccupiedOrbital) +
                                                     2 * δ(p, s) * δ(q, r) * constrain(p => OccupiedOrbital, q => VirtualOrbital)

    # Issue #13 from ExcitationOperators.jl
    e1 = hf_expectation_value(∑(E(k, j) * E(j, i), [i, j, k]))
    e2 = hf_expectation_value(∑(E(i, j) * E(j, i), [i, j]))
    @test simplify(e1) == simplify(e2)

    @test hf_expectation_value(E(i, j) * E(j, i) * constrain(i => OccupiedOrbital, j => OccupiedOrbital)) == 4 * δ(i, j) * constrain(i => OccupiedOrbital)
end

@testset "jacobi_identity" begin
    p = general(1)
    q = general(2)
    r = general(3)
    s = general(4)
    t = general(5)
    u = general(6)

    f = commutator
    A = E(p, q)
    B = E(r, s)
    C = E(t, u)

    @test f(A, f(B, C)) + f(C, f(A, B)) + f(B, f(C, A)) == SASQ.Expression(0)
end

@testset "psym_tensor" begin
    p = general(1)
    q = general(2)
    r = general(3)
    s = general(4)
    i = occupied(1)
    j = occupied(2)
    k = occupied(3)
    a = virtual(1)
    b = virtual(2)
    c = virtual(3)

    @test psym_tensor("g", p, q, r, s) == psym_tensor("g", r, s, p, q)
    @test psym_tensor("g", a, i, b, j) == psym_tensor("g", b, j, a, i)
    @test psym_tensor("g", p, q, b, j) == psym_tensor("g", b, j, p, q)
    @test psym_tensor("t", a, i, b, j, c, k) == psym_tensor("t", b, j, c, k, a, i)
end

@testset "hf_equations" begin
    p = general(1)
    q = general(2)
    r = general(3)
    s = general(4)
    i = occupied(1)
    j = occupied(2)

    h = ∑(rsym_tensor("h", p, q) * E(p, q), [p, q])
    g = 1 // 2 * ∑(rsym_tensor("g", p, q, r, s) * e(p, q, r, s), [p, q, r, s])
    H = simplify(h + g)

    energy = simplify(hf_expectation_value(H))
    energy_exp = simplify(2 // 1 * ∑(rsym_tensor("h", p, p) * constrain(p => OccupiedOrbital), [p]) + ∑((2 // 1 * rsym_tensor("g", p, p, q, q) - rsym_tensor("g", p, q, q, p)) * constrain(p => OccupiedOrbital, q => OccupiedOrbital), [p, q]))
    @test energy == energy_exp

    gradient = simplify(hf_expectation_value(commutator(H, E(p, q) - E(q, p))))
    fock = (rsym_tensor("h", p, q) + ∑(constrain(i => OccupiedOrbital) * (2 // 1 * rsym_tensor("g", p, q, i, i) - rsym_tensor("g", p, i, i, q)), [i])) * constrain(p => VirtualOrbital, q => OccupiedOrbital)

    @test gradient * constrain(p => OccupiedOrbital, q => OccupiedOrbital) == SASQ.Expression(0)
    @test gradient * constrain(p => VirtualOrbital, q => VirtualOrbital) == SASQ.Expression(0)
    @test simplify(gradient * constrain(p => VirtualOrbital, q => OccupiedOrbital)) == simplify(4 // 1 * fock * constrain(p => VirtualOrbital, q => OccupiedOrbital))
end


@testset "ccsd_equations" begin
    p = general(1)
    q = general(2)
    r = general(3)
    s = general(4)

    h = ∑(rsym_tensor("h", p, q) * E(p, q), [p, q])
    g = 1 // 2 * ∑(rsym_tensor("g", p, q, r, s) * e(p, q, r, s), [p, q, r, s])
    H = simplify(h + g)

    T1 = ∑(psym_tensor("t", p, q) * constrain(p => VirtualOrbital, q => OccupiedOrbital) * E(p, q), [p, q])
    T2 = ∑(psym_tensor("t", p, q, r, s) * constrain(p => VirtualOrbital, q => OccupiedOrbital, r => VirtualOrbital, s => OccupiedOrbital) * E(p, q) * E(r, s), [p, q, r, s])

    #@show simplify(hf_expectation_value(H * T1)) == SASQ.Expression(0)
end