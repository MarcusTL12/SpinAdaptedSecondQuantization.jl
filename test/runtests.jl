using Test

using SpinAdaptedSecondQuantization

@testset "indices" begin
    println()

    p = gen("p")
    i = occ("i")
    a = vir("a")
    @show p i a

    @test (space(p) <: GeneralOrbital) == true
    @test (space(p) <: OccupiedOrbital) == false
    @test (space(i) <: GeneralOrbital) == true
    @test (space(i) <: OccupiedOrbital) == true
    
    @test (space(i) == GeneralOrbital) == false
    @test (space(i) == OccupiedOrbital) == true

    @test isdisjoint(p, i) == false
    @test isdisjoint(a, i) == true

    @test isocc(i) == true
    @test isvir(i) == false
    @test isocc(a) == false
    @test isvir(a) == true

    println()
end

@testset "kroenecker delta" begin
    println()

    p = gen("p")
    i = occ("i")
    a = vir("a")

    dpa = SpinAdaptedSecondQuantization.KroeneckerDelta(p, a)
    dia = SpinAdaptedSecondQuantization.KroeneckerDelta(i, a)
    dip = SpinAdaptedSecondQuantization.KroeneckerDelta(i, p)

    @show dpa dia dip
    
    println()
end
