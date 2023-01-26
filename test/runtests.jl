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

    @test string(a) == "\e[36ma\e[39m"

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

@testset "singlet excitation operator" begin
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

@testset "real tensor" begin
    println()

    p = general(1)
    q = general(2)

    hpq = SpinAdaptedSecondQuantization.RealTensor("h", [p, q])

    @show hpq

    @test string(hpq) == "h_pq"

    println()
end

@testset "term show" begin
    println()

    p = general(1)
    q = general(2)
    i = occupied(1)
    a = virtual(1)

    t = SpinAdaptedSecondQuantization.Term(
        3 // 5,
        [i, a],
        [SpinAdaptedSecondQuantization.KroneckerDelta(p, a)],
        [
            SpinAdaptedSecondQuantization.RealTensor("h", [p, a]),
            SpinAdaptedSecondQuantization.RealTensor("g", [i, a, i, i])
        ],
        [SpinAdaptedSecondQuantization.SingletExcitationOperator(p, q)],
    )

    @show t
    @show SpinAdaptedSecondQuantization.get_all_indices(t)

    println()
end

@testset "expression constructors" begin
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
    @test string(dai) == "0"

    println()
end

@testset "term exchange_indices" begin
    println()

    p = general(1)
    q = general(2)
    i = occupied(1)
    j = occupied(2)
    a = virtual(1)
    b = virtual(2)

    t = SpinAdaptedSecondQuantization.Term(
        3 // 5,
        [i, a],
        [SpinAdaptedSecondQuantization.KroneckerDelta(p, a)],
        [
            SpinAdaptedSecondQuantization.RealTensor("h", [p, a]),
            SpinAdaptedSecondQuantization.RealTensor("g", [i, a, i, q])
        ],
        [SpinAdaptedSecondQuantization.SingletExcitationOperator(p, q)],
    )

    @show t

    t2 = SpinAdaptedSecondQuantization.exchange_indices(
        t,
        [p => q, q => p, a => b, i => j]
    )

    @test string(t2) == "3/5 ∑_\e[92mj\e[39m\e[36mb\e[39m(δ_q\e[36mb\e[39m \
g_\e[92mj\e[39m\e[36mb\e[39m\e[92mj\e[39mp h_q\e[36mb\e[39m E_qp)"

    @show t2

    t3 = SpinAdaptedSecondQuantization.make_space_for_indices(t, i)
    @show t3

    t4 = SpinAdaptedSecondQuantization.make_space_for_indices(t, [i, j, a])
    @show t4

    println()
end

@testset "term summation" begin
    println()

    p = general(1)
    q = general(2)
    i = occupied(1)
    j = occupied(2)
    a = virtual(1)
    b = virtual(2)

    t = SpinAdaptedSecondQuantization.Term(
        1,
        SpinAdaptedSecondQuantization.MOIndex[],
        [SpinAdaptedSecondQuantization.KroneckerDelta(p, a)],
        [SpinAdaptedSecondQuantization.RealTensor("h", [p, q])],
        SpinAdaptedSecondQuantization.Operator[],
    )

    @show t

    t = ∑(t, [p, q])

    @show t

    println()
end

@testset "term summation simplify" begin
    println()

    p = general(1)
    q = general(2)
    i = occupied(1)
    j = occupied(2)
    a = virtual(1)
    b = virtual(2)

    t = real_tensor("g", b, i, a, j)

    @show t

    t = ∑(t, [a, b, i, j])

    @show t

    @show SpinAdaptedSecondQuantization.sort_summation_indices(t.terms[1])

    println()
end
