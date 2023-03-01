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
