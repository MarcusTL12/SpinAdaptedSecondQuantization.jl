using Test

using Base.Threads: @spawn

using SpinAdaptedSecondQuantization

tasks = Task[]

macro testmacro(code)
    push!(tasks, @spawn begin
        tmppath, io = mktemp()
        redirect_stdout(io) do
            eval(code)
        end
        close(io)
        read(tmppath, String)
    end)
end

@testmacro @testset "indices" begin
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

    @test string(a) != "a"

    println()
end

@testmacro @testset "kronecker delta" begin
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

@testmacro @testset "singlet excitation operator" begin
    println()

    p = general(1)
    q = general(2)
    i = occupied(1)
    a = virtual(1)

    Epq = SpinAdaptedSecondQuantization.SingletExcitationOperator(p, q)
    Eai = SpinAdaptedSecondQuantization.SingletExcitationOperator(a, i)
    Eip = SpinAdaptedSecondQuantization.SingletExcitationOperator(i, p)

    @show Epq Eai Eip

    @test string(Epq) == "E_pq"

    println()
end

@testmacro @testset "real tensor" begin
    println()

    p = general(1)
    q = general(2)

    hpq = SpinAdaptedSecondQuantization.RealTensor("h", [p, q])

    @show hpq

    @test string(hpq) == "h_pq"

    println()
end

@testmacro @testset "term" begin
    println()

    p = general(1)
    q = general(2)
    i = occupied(1)
    a = virtual(1)

    Epq = E(p, q)
    dpi = δ(p, i)
    dai = δ(a, i)
    hai = real_tensor("h", a, i)

    @show Epq dpi dai hai

    @test string(Epq) == "E_pq"

    println()
end

for t in tasks
    println(fetch(t))
end
