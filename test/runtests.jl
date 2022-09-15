using Test

using SpinAdaptedSecondQuantization

@testset "indices" begin
    p = gen("p")
    i = occ("i")
    a = vir("a")
    @show p i a

    @test (p.space <: GeneralOrbital) == true
    @test (p.space <: OccupiedOrbital) == false
    @test (i.space <: GeneralOrbital) == true
    @test (i.space <: OccupiedOrbital) == true
    
    @test (i.space == GeneralOrbital) == false
    @test (i.space == OccupiedOrbital) == true

    @test isdisjoint(p, i) == false
    @test isdisjoint(a, i) == true
end
