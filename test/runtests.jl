using Test

using SpinAdaptedSecondQuantization

@testset "greet" begin
    @test 1 == 1
    @test 2 == 1 broken=true
    SpinAdaptedSecondQuantization.greet()
end
