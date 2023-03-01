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

    @test string(t) == "3/5 ∑_rs(δ_pqs g_rsrr h_ps E_pq) C(p∈V, r∈O)"

    t = SASQ.lower_delta_indices(t)
    @test string(t) == "3/5 ∑_rs(δ_pqs g_rprr h_pp E_pp) C(p∈V, r∈O)"
end
