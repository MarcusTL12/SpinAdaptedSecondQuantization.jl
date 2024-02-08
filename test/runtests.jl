using Test

using SpinAdaptedSecondQuantization
SASQ.disable_color()

@testset "ortital spaces" begin
    @test OccupiedOrbital <: GeneralOrbital
    @test VirtualOrbital <: GeneralOrbital
    @test isdisjoint(VirtualOrbital, OccupiedOrbital)
    @test !isdisjoint(GeneralOrbital, OccupiedOrbital)
    @test typeintersect(GeneralOrbital, OccupiedOrbital) == OccupiedOrbital
    @test typeintersect(OccupiedOrbital, VirtualOrbital) == Union{}
    @test SASQ.is_strict_subspace(OccupiedOrbital, GeneralOrbital)
    @test !SASQ.is_strict_subspace(GeneralOrbital, GeneralOrbital)
end

@testset "kronecker delta" begin
    dpp = δ(1, 1)
    dpq = δ(1, 2)

    @test string(dpp) == "1"
    @test string(dpq) == "δ_pq"
    @test dpp == SASQ.Expression(1)
end

@testset "kronecker delta compact" begin
    dpq = SASQ.KroneckerDelta(1, 2)
    dqr = SASQ.KroneckerDelta(2, 3)
    dpqr = SASQ.KroneckerDelta(1, 2, 3)

    @test SASQ.compact_deltas([dpq, dqr]) == [dpqr]
end

@testset "print singlet excitation operator" begin
    Epq = E(1, 2)
    @test string(Epq) == "E_pq"
end

@testset "print real tensor" begin
    hpq = real_tensor("h", 1, 2)
    @test string(hpq) == "h_pq"

    gpqrs = real_tensor("g", 1, 2, 3, 4)
    @test string(gpqrs) == "g_pqrs"
end

@testset "print term" begin
    p = 1
    q = 2
    i = 3
    a = 4

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

    @test string(SASQ.Expression([t])) ==
          "3/5 ∑_ia(δ_pqa g_iaii h_pa E_pq) C(p∈v, q∈v)"

    t = SASQ.lower_delta_indices(t)
    @test string(SASQ.Expression([t])) ==
          "3/5 ∑_ia(δ_pqa g_ipii h_pp E_pp) C(p∈v, q∈v)"

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

    @test string(SASQ.Expression([t])) ==
          "3/5 ∑_ia(δ_pa g_iaiq h_pa E_pq) C(p∈v)"
end

@testset "term exchange_indices" begin
    p = 1
    q = 2
    r = 3
    i = 4
    a = 5

    t = SASQ.Term(
        3 // 7,
        Int[],
        [SASQ.KroneckerDelta(p, q, r)],
        [
            SASQ.RealTensor("h", [p, a]),
            SASQ.RealTensor("g", [i, a, i, q])
        ],
        [SASQ.SingletExcitationOperator(p, q)],
        SASQ.Constraints(i => OccupiedOrbital, a => VirtualOrbital)
    )

    @test string(SASQ.Expression([t])) ==
          "3/7 δ_pqr g_stsq h_pt E_pq C(s∈o, t∈v)"

    t = SASQ.lower_delta_indices(t)

    @test string(SASQ.Expression([t])) ==
          "3/7 δ_pqr g_stsp h_pt E_pp C(s∈o, t∈v)"

    t2 = SASQ.exchange_indices(t, [p => q])

    @test string(SASQ.Expression([t2])) ==
          "3/7 δ_qr g_stsq h_qt E_qq C(s∈o, t∈v)"

    t3 = SASQ.exchange_indices(t2, [q => r])

    @test string(SASQ.Expression([t3])) == "3/7 g_stsr h_rt E_rr C(s∈o, t∈v)"
end

@testset "term summation delta" begin
    p = 1
    q = 2
    i = 3
    a = 4

    t = SASQ.Term(
        1,
        Int[],
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
    p = 1
    q = 2
    r = 3
    s = 4

    t1 = SASQ.Term(
        3 // 5,
        [p, q],
        SASQ.KroneckerDelta[],
        SASQ.Tensor[],
        [SASQ.SingletExcitationOperator(p, q)],
        SASQ.Constraints(p => OccupiedOrbital, q => VirtualOrbital)
    )

    @test simplify(t1 * t1) == SASQ.Term(9 // 25,
        [p, q, r, s],
        SASQ.KroneckerDelta[],
        SASQ.Tensor[],
        [
            SASQ.SingletExcitationOperator(p, q),
            SASQ.SingletExcitationOperator(r, s)
        ],
        SASQ.Constraints(
            p => OccupiedOrbital,
            q => VirtualOrbital,
            r => OccupiedOrbital,
            s => VirtualOrbital
        )
    )
end

@testset "expression addition" begin
    p = 1
    q = 2
    r = 3
    s = 4

    i = 5
    j = 6
    ic = constrain(i => OccupiedOrbital)
    jc = constrain(j => OccupiedOrbital)

    a = 7
    b = 8
    ac = constrain(a => VirtualOrbital)
    bc = constrain(b => VirtualOrbital)

    e1 = δ(a, b) * ac * bc
    e2 = real_tensor("h", i, j) * ic * jc

    @test (e1 + 1) == (1 + e1)
    @test (e1 + e2) == (e2 + e1)
    @test (e1 + e2 + 2) == (2 + e1 + e2) == (e2 + 2 + e1)
end

@testset "expression multiplication" begin
    p = 1
    q = 2
    r = 3
    s = 4

    i = 5
    j = 6
    k = 7
    l = 8
    ic = constrain(i => OccupiedOrbital)
    jc = constrain(j => OccupiedOrbital)
    kc = constrain(k => OccupiedOrbital)
    lc = constrain(l => OccupiedOrbital)

    m = 9
    n = 10

    a = 11
    b = 12
    ac = constrain(a => VirtualOrbital)
    bc = constrain(b => VirtualOrbital)

    e1 = δ(a, b) * ac * bc
    e2 = real_tensor("h", i, j) * ic * jc
    h = ∑(e2 * E(i, j), [i, j])
    g = ∑(
        real_tensor("g", i, j, k, l) * E(i, j) * E(k, l) * ic * jc * kc * lc,
        [i, j, k, l]
    )

    @test ((e1 + e2) * 3) // 5 * e2 ==
          3 // 5 * δ(a, b) * real_tensor("h", i, j) * ac * bc * ic * jc +
          3 // 5 * real_tensor("h", i, j) * real_tensor("h", i, j) * ic * jc

    @test (h + h) == 2 * h

    @test simplify(g * h) == ∑(
        real_tensor("g", 1, 2, 3, 4) *
        real_tensor("h", 5, 6) * E(1, 2) * E(3, 4) * E(5, 6) *
        occupied(1:6...),
        collect(1:6)
    )

    @test simplify(h * g) == ∑(
        real_tensor("g", 3, 4, 5, 6) *
        real_tensor("h", 1, 2) * E(1, 2) * E(3, 4) * E(5, 6) *
        occupied(1:6...),
        collect(1:6)
    )
end

@testset "simple commutator" begin
    p = 1
    q = 2
    r = 3
    s = 4

    i = 5
    j = 6
    k = 7

    a = 9
    b = 10

    Eai = E(a, i) * occupied(i) * virtual(a)
    Eia = E(i, a) * occupied(i) * virtual(a)

    @test commutator(E(p, q), E(r, s)) == E(p, s) * δ(q, r) - E(r, q) * δ(p, s)
    @test commutator(Eia, Eai) == (E(i, i) - E(a, a)) * occupied(i) * virtual(a)
    @test commutator(E(a, i), E(a, i)) == SASQ.Expression(0)

    h = ∑(real_tensor("h", p, q) * E(p, q), [p, q])
    hcomm = simplify(commutator(h, E(p, q)))
    @test hcomm == ∑(
        -real_tensor("h", q, r) * E(p, r) +
        real_tensor("h", r, p) * E(r, q),
        [r]
    )
    @test hcomm * E(r, r) == ∑(
        -real_tensor("h", q, s) * E(p, s) * E(r, r) +
        real_tensor("h", s, p) * E(s, q) * E(r, r),
        [s]
    )
end

@testset "hf_expectation_value" begin
    @test hf_expectation_value(E(1, 2)) == 2δ(1, 2) * occupied(1, 2)
    @test hf_expectation_value(E(1, 2) * virtual(1)) |> iszero

    # Issue #13 from ExcitationOperators.jl
    e1 = hf_expectation_value(∑(E(3, 2) * E(2, 1), 1:3))
    e2 = hf_expectation_value(∑(E(1, 2) * E(2, 1), 1:2))
    @test simplify(e1) == simplify(e2)

    @test hf_expectation_value(E(1, 2) * E(2, 1) * occupied(1, 2)) ==
          4δ(1, 2) * occupied(1, 2)
end

@testset "jacobi_identity" begin
    f = commutator
    A = E(1, 2)
    B = E(3, 4)
    C = E(5, 6)

    @test iszero(f(A, f(B, C)) + f(C, f(A, B)) + f(B, f(C, A)))
end

@testset "hf energy" begin
    h = ∑(real_tensor("h", 1, 2) * E(1, 2), 1:2)
    g = 1 // 2 * ∑(real_tensor("g", 1:4...) * e(1:4...), 1:4) |> simplify

    H = simplify(h + g)

    E_hf = H |> hf_expectation_value |> simplify
    @test E_hf == ∑(2real_tensor("h", 1, 1) * occupied(1), [1]) + ∑(
        (2real_tensor("g", 1, 1, 2, 2) - real_tensor("g", 1, 2, 2, 1)) *
        occupied(1, 2),
        1:2
    )

    hF = ∑(
        (real_tensor("F", 1, 2) + ∑(
            (-2real_tensor("g", 1, 2, 3, 3) + real_tensor("g", 1, 3, 3, 2)) *
            occupied(3),
            [3]
        )) * E(1, 2),
        1:2
    )
    HF = hF + g
    E_hf = ((HF + H) // 2) |> hf_expectation_value |> simplify
    @test E_hf == ∑(
        (real_tensor("h", 1, 1) + real_tensor("F", 1, 1)) * occupied(1),
        [1]
    )
end

@testset "wicks theorem" begin
    @test SASQ.fully_contracted_pairs([1, 2]) == [[(1, 2)]]
    @test SASQ.fully_contracted_pairs([1, 2, 3, 4]) == [[(1, 2), (3, 4)],
        [(1, 3), (2, 4)],
        [(1, 4), (2, 3)]]
    @test SASQ.fully_contracted_pairs([1, 2, 3, 4, 5, 6]) ==
          [
        [(1, 2), (3, 4), (5, 6)],
        [(1, 2), (3, 5), (4, 6)],
        [(1, 2), (3, 6), (4, 5)],
        [(1, 3), (2, 4), (5, 6)],
        [(1, 3), (2, 5), (4, 6)],
        [(1, 3), (2, 6), (4, 5)],
        [(1, 4), (2, 3), (5, 6)],
        [(1, 4), (2, 5), (3, 6)],
        [(1, 4), (2, 6), (3, 5)],
        [(1, 5), (2, 3), (4, 6)],
        [(1, 5), (2, 4), (3, 6)],
        [(1, 5), (2, 6), (3, 4)],
        [(1, 6), (2, 3), (4, 5)],
        [(1, 6), (2, 4), (3, 5)],
        [(1, 6), (2, 5), (3, 4)]
    ]

    h = sum(
        ∑(real_tensor("h", 1, 2) * fermiondag(1, spin) * fermion(2, spin), 1:2)
        for spin in (α, β)
    )
    g = 1 // 2 * sum(
        ∑(real_tensor("g", 1:4...) *
          fermiondag(1, spin1) *
          fermiondag(3, spin2) *
          fermion(4, spin2) *
          fermion(2, spin1), 1:4)
        for spin1 in (α, β), spin2 in (α, β)
    )
    H = h + g
    E_hf = simplify(wick_theorem(H))
    @test E_hf == ∑(2 * real_tensor("h", 1, 1) * occupied(1), [1]) +
                  ∑((2 * real_tensor("g", 1, 1, 2, 2) -
                     real_tensor("g", 1, 2, 2, 1)) * occupied(1, 2), 1:2)
end

@testset "cis" begin
    hF = ∑(
        (real_tensor("F", 1, 2) + ∑(
            (-2psym_tensor("g", 1, 2, 3, 3) + psym_tensor("g", 1, 3, 3, 2)) *
            occupied(3),
            [3]
        )) * E(1, 2),
        1:2
    )
    g = 1 // 2 * ∑(psym_tensor("g", 1:4...) * e(1:4...), 1:4) |> simplify

    HF = simplify(hF + g)

    EHF = hf_expectation_value(HF) |> simplify

    E_single = simplify(
        (hf_expectation_value(E(1, 2) * HF * E(2, 1)) // 2 -
         hf_expectation_value(HF)) * occupied(1) * virtual(2)
    )

    @test E_single ==
          (real_tensor("F", 2, 2) - real_tensor("F", 1, 1) +
           2psym_tensor("g", 1, 2, 2, 1) - psym_tensor("g", 1, 1, 2, 2)) *
          occupied(1) * virtual(2)

    H_cis = simplify_heavy(
        (hf_expectation_value(E(2, 3) * HF * E(4, 1)) // 2 -
         hf_expectation_value(HF) * δ(1, 2) * δ(3, 4)) *
        occupied(1, 2) * virtual(3, 4)
    )

    @test H_cis ==
          (δ(1, 2) * real_tensor("F", 3, 4) - δ(3, 4) * real_tensor("F", 1, 2) +
           2psym_tensor("g", 1, 4, 3, 2) - psym_tensor("g", 1, 2, 3, 4)) *
          occupied(1, 2) * virtual(3, 4)
end

@testset "simplify_heavy" begin
    x = simplify(∑(real_tensor("g", 1, 2) * real_tensor("g", 2, 1) *
                   (real_tensor("h", 1, 2) - real_tensor("h", 2, 1)), 1:2))

    @test !iszero(x)
    @test iszero(simplify_heavy(x))
end

@testset "ket" begin
    h = ∑(real_tensor("h", 1, 2) * E(1, 2), 1:2)
    g = 1 // 2 * ∑(real_tensor("g", 1:4...) * e(1:4...), 1:4)
    H = simplify(h + g)

    Hket = act_on_ket(H) |> simplify
    noop_terms = filter(x -> iszero(length(x.operators)), Hket.terms)
    Hket0 = SASQ.Expression(noop_terms)

    @test Hket0 == ∑(2real_tensor("h", 1, 1) * occupied(1), [1]) + ∑(
        (2real_tensor("g", 1, 1, 2, 2) - real_tensor("g", 1, 2, 2, 1)) *
        occupied(1, 2),
        1:2
    )
end

@testset "boson" begin
    b = boson()
    bdag = bosondag()
    a = fermion(1, α)
    adag = fermiondag(1, α)

    @test act_on_ket(b * bdag) == SASQ.Expression(1)
    @test act_on_ket(b * b * bdag * bdag) == SASQ.Expression(2)
    @test act_on_ket(b * b * b * bdag * bdag * bdag) == SASQ.Expression(6)
    @test act_on_ket(a * b * adag * bdag) == virtual(1)
    @test simplify(a * b * bdag * adag * E(2, 3)) == a * adag * E(2, 3) * b * bdag
end

@testset "PermuteTensor" begin
    @test permute_tensor(1, 2, 3, 4) == permute_tensor(3, 4, 1, 2)
    eq1 = permute_tensor(1, 2, 3, 4) * real_tensor("h", 3, 4)
    eq2 = permute_tensor(1, 2, 3, 4) * real_tensor("h", 1, 2)
    @test simplify_heavy(eq1) == eq2
end

@testset "print code" begin
    equation = summation(real_tensor("h", 1, 2) * real_tensor("g", 2, 1, 3) * occupied(1) * virtual(2,3), 1:2)
    trans = translate(VirtualOrbital => [3])

    code = print_code(equation.terms[1], "omega", trans)
    expected_code = """omega_a +=  +1.00000000 * np.einsum("bia,ib->a", g_vov, h_ov, optimize="optimal");"""
    @test code == expected_code

    code_eT = print_eT_code(equation.terms[1], "omega", trans, "test")
    expected_code_eT = """print(generate_eT_code_from_einsum(\n    routine_name=\"test\",\n    prefactor= +1.00000000,\n    contraction_string=\"bia,ib->a\",\n    arrays=[g_vov, h_ov, omega],\n    symbols=[\"g_vov\", \"h_ov\", \"omega\"],\n), end='!\\n!\\n')"""
    @test code_eT == expected_code_eT
end
