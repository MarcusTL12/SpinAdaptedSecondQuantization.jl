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
    q = general(2)
    i = occupied(1)
    a = virtual(1)

    dpa = SpinAdaptedSecondQuantization.KroneckerDelta(p, a)
    dia = SpinAdaptedSecondQuantization.KroneckerDelta(i, a)
    dip = SpinAdaptedSecondQuantization.KroneckerDelta(i, p)
    dpp = SpinAdaptedSecondQuantization.KroneckerDelta(p, p)
    dpqi = SpinAdaptedSecondQuantization.KroneckerDelta(p, q, i)

    @show dpa dia dip dpp dpqi

    @test dpa != 0
    @test dia == 0
    @test dip != 0
    @test dpp == 1
    @test dpqi != 0

    println()
end

@testset "kronecker delta compact" begin
    println()

    p = general(1)
    q = general(2)
    r = general(3)
    i = occupied(1)
    a = virtual(1)

    dpa = SpinAdaptedSecondQuantization.KroneckerDelta(p, a)
    dri = SpinAdaptedSecondQuantization.KroneckerDelta(r, i)
    dpq = SpinAdaptedSecondQuantization.KroneckerDelta(p, q)
    dqr = SpinAdaptedSecondQuantization.KroneckerDelta(q, r)

    v = [dpa, dpq, dri]

    @show v

    v = SpinAdaptedSecondQuantization.compact_deltas(v)

    @show v
    @test v != 0

    v = [dpa, dpq, dri, dqr]

    @show v

    v = SpinAdaptedSecondQuantization.compact_deltas(v)

    @show v
    @test v == 0

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
        [SpinAdaptedSecondQuantization.KroneckerDelta(p, a),
            SpinAdaptedSecondQuantization.KroneckerDelta(a, q)],
        [
            SpinAdaptedSecondQuantization.RealTensor("h", [p, a]),
            SpinAdaptedSecondQuantization.RealTensor("g", [i, a, i, i])
        ],
        [SpinAdaptedSecondQuantization.SingletExcitationOperator(p, q)]
    )

    @show t
    @show SpinAdaptedSecondQuantization.get_all_indices(t)

    t = SpinAdaptedSecondQuantization.lower_delta_indices(t)

    @show t

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
    r = general(3)
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

    t = SpinAdaptedSecondQuantization.Term(
        3 // 7,
        SpinAdaptedSecondQuantization.MOIndex[],
        [SpinAdaptedSecondQuantization.KroneckerDelta(p, q, r)],
        [
            SpinAdaptedSecondQuantization.RealTensor("h", [p, a]),
            SpinAdaptedSecondQuantization.RealTensor("g", [i, a, i, q])
        ],
        [SpinAdaptedSecondQuantization.SingletExcitationOperator(p, q)],
    )

    println()

    @show t

    t = SpinAdaptedSecondQuantization.lower_delta_indices(t)

    @show t

    t2 = SpinAdaptedSecondQuantization.exchange_indices(t, [p => q])

    @show t2

    t3 = SpinAdaptedSecondQuantization.exchange_indices(t2, [q => r])

    @show t3

    t4 = SpinAdaptedSecondQuantization.exchange_indices(t, [q => p])

    println()
    @show t4

    t5 = SpinAdaptedSecondQuantization.exchange_indices(t4, [r => p])

    @show t5

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
        [SpinAdaptedSecondQuantization.KroneckerDelta(a, p)],
        [SpinAdaptedSecondQuantization.RealTensor("h", [a, i])],
        SpinAdaptedSecondQuantization.Operator[
            SpinAdaptedSecondQuantization.SingletExcitationOperator(p, q)
        ],
    )

    @show t

    t = SpinAdaptedSecondQuantization.lower_delta_indices(t)

    @show t

    println("\np -> a:")
    t1 = ∑(t, [p])
    @show t1
    t1 = SpinAdaptedSecondQuantization.simplify_summation_deltas(t1)
    @show t1
    t2 = ∑(t1, [a])
    @show t2

    println("\na -> p:")
    t1 = ∑(t, [a])
    @show t1
    t1 = SpinAdaptedSecondQuantization.simplify_summation_deltas(t1)
    @show t1
    t2 = ∑(t1, [p])
    @show t2

    println("\nap:")
    t1 = ∑(t, [a, p])
    @show t1
    t1 = SpinAdaptedSecondQuantization.simplify_summation_deltas(t1)
    @show t1

    println()
end

@testset "term summation simplify" begin
    println()

    p = general(1)
    q = general(2)
    i = occupied(1)
    j = occupied(2)
    b = virtual(2)
    d = virtual(4)

    t = real_tensor("g", b, i, d, j)

    @show t

    t = ∑(t, [b, d, i, j])

    @show t

    t = SpinAdaptedSecondQuantization.lower_summation_indices(t.terms[1])

    @show t

    t = SpinAdaptedSecondQuantization.sort_summation_indices(t)

    @show t

    println()
end

@testset "term exchange_indices constraints" begin
    println()

    p = general(1)
    q = general(2)
    r = general(3)
    s = general(4)
    i = occupied(1)
    j = occupied(2)
    a = virtual(1)
    b = virtual(2)

    t = SpinAdaptedSecondQuantization.Term(
        3 // 5,
        SpinAdaptedSecondQuantization.MOIndex[],
        SpinAdaptedSecondQuantization.KroneckerDelta[],
        [
            SpinAdaptedSecondQuantization.RealTensor("h", [a, b]),
        ],
        [SpinAdaptedSecondQuantization.SingletExcitationOperator(i, j)],
    )

    @show t

    t2 = SpinAdaptedSecondQuantization.exchange_indices(
        t,
        [a => p, b => q, i => r, j => s]
    )

    @show t2

    t3 = SpinAdaptedSecondQuantization.exchange_indices(
        t2,
        [p => b]
    )

    @show t3

    t4 = ∑(t3, [s])

    @show t4

    println()
end
