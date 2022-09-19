using Test

using SpinAdaptedSecondQuantization

@testset "indices" begin
    println()

    p = general(1)
    i = occupied(1)
    a = virtual(1)
    p1 = general(8)
    @show p i a p1

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

    println()
end

@testset "kronecker delta" begin
    println()

    p = general(1)
    i = occupied(1)
    a = virtual(1)

    dpa = SpinAdaptedSecondQuantization.KroneckerDelta(p, a)
    dia = SpinAdaptedSecondQuantization.KroneckerDelta(i, a)
    dip = SpinAdaptedSecondQuantization.KroneckerDelta(i, p)

    @show dpa dia dip

    @test dpa != 0
    @test dia == 0
    @test dip != 0

    println()
end
