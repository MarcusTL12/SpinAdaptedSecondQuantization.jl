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
    @test SASQ.is_strict_subspace(OccupiedOrbital, GeneralOrbital)
    @test !SASQ.is_strict_subspace(GeneralOrbital, GeneralOrbital)
end

@testset "print indices" begin
    @test SASQ.print_mo_index(1) == "p"
    @test SASQ.print_mo_index(2) == "q"
    @test SASQ.print_mo_index(9) == "p₁"
    @test SASQ.print_mo_index(104) == "w₁₂"
end

@testset "kronecker delta" begin
    dpp = SASQ.KroneckerDelta(1, 1)
    dpq = SASQ.KroneckerDelta(1, 2)

    @test string(dpp) == "1"
    @test string(dpq) == "δ_pq"
    @test dpp == 1
end

@testset "kronecker delta compact" begin
    dpq = SASQ.KroneckerDelta(1, 2)
    dqr = SASQ.KroneckerDelta(2, 3)
    dpqr = SASQ.KroneckerDelta(1, 2, 3)

    @test SASQ.compact_deltas([dpq, dqr]) == [dpqr]
end

@testset "print singlet excitation operator" begin
    Epq = SASQ.SingletExcitationOperator(1, 2)
    @test string(Epq) == "E_pq"
end

@testset "print real tensor" begin
    hpq = SASQ.RealTensor("h", [1, 2])
    @test string(hpq) == "h_pq"

    gpqrs = SASQ.RealTensor("g", [1, 2, 3, 4])
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

    @test string(t) == "3/5 ∑_rs(δ_pqs g_rsrr h_ps E_pq) C(p∈V, q∈V, r∈O, s∈V)"

    t = SASQ.lower_delta_indices(t)
    @test string(t) == "3/5 ∑_rs(δ_pqs g_rprr h_pp E_pp) C(p∈V, q∈V, r∈O, s∈V)"

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

    @test string(t) == "3/5 ∑_rs(δ_ps g_rsrq h_ps E_pq) C(p∈V, r∈O, s∈V)"
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

    @test string(t) == "3/7 δ_pqr g_stsq h_pt E_pq C(s∈O, t∈V)"

    t = SASQ.lower_delta_indices(t)

    @test string(t) == "3/7 δ_pqr g_stsp h_pt E_pp C(s∈O, t∈V)"

    t2 = SASQ.exchange_indices(t, [p => q])

    @test string(t2) == "3/7 δ_qr g_stsq h_qt E_qq C(s∈O, t∈V)"

    t3 = SASQ.exchange_indices(t2, [q => r])

    @test string(t3) == "3/7 g_stsr h_rt E_rr C(s∈O, t∈V)"
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
        SASQ.Constraints(i=>OccupiedOrbital, a => VirtualOrbital)
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
