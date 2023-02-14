using Test

using SpinAdaptedSecondQuantization

@testset "indices" begin
    p = general(1)
    i = occupied(1)
    a = virtual(1)
    p1 = general(8)

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
    @test string(i) == "\e[92mi\e[39m"
end

@testset "kronecker delta" begin
    p = general(1)
    q = general(2)
    i = occupied(1)
    a = virtual(1)

    dpa = SpinAdaptedSecondQuantization.KroneckerDelta(p, a)
    dia = SpinAdaptedSecondQuantization.KroneckerDelta(i, a)
    dip = SpinAdaptedSecondQuantization.KroneckerDelta(i, p)
    dpp = SpinAdaptedSecondQuantization.KroneckerDelta(p, p)
    dpqi = SpinAdaptedSecondQuantization.KroneckerDelta(p, q, i)

    @test string(dpa) == "δ_p\e[36ma\e[39m"
    @test string(dia) == "0"
    @test string(dip) == "δ_p\e[92mi\e[39m"
    @test string(dpp) == "1"
    @test string(dpqi) == "δ_pq\e[92mi\e[39m"

    @test dpa != 0
    @test dia == 0
    @test dip != 0
    @test dpp == 1
    @test dpqi != 0
end

@testset "kronecker delta compact" begin
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

    @test string(v) == "KroneckerDelta[δ_p\e[36ma\e[39m, δ_pq, δ_r\e[92mi\e[39m]"

    v = SpinAdaptedSecondQuantization.compact_deltas(v)

    @test string(v) == "KroneckerDelta[δ_pq\e[36ma\e[39m, δ_r\e[92mi\e[39m]"
    @test v != 0

    v = [dpa, dpq, dri, dqr]

    @test string(v) == "KroneckerDelta[δ_p\e[36ma\e[39m, δ_pq, δ_r\e[92mi\e[39m, δ_qr]"

    v = SpinAdaptedSecondQuantization.compact_deltas(v)

    @test string(v) == "0"
    @test v == 0
end

@testset "singlet excitation operator" begin
    p = general(1)
    q = general(2)
    i = occupied(1)
    a = virtual(1)

    Epq = SpinAdaptedSecondQuantization.SingletExcitationOperator(p, q)
    Eai = SpinAdaptedSecondQuantization.SingletExcitationOperator(a, i)
    Eip = SpinAdaptedSecondQuantization.SingletExcitationOperator(i, p)

    @test string(Epq) == "E_pq"
    @test string(Eai) == "E_\e[36ma\e[39m\e[92mi\e[39m"
    @test string(Eip) == "E_\e[92mi\e[39mp"
end

@testset "real tensor" begin
    p = general(1)
    q = general(2)

    hpq = SpinAdaptedSecondQuantization.RealTensor("h", [p, q])

    @test string(hpq) == "h_pq"
end

@testset "term show" begin
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

    @test string(t) == "3/5 ∑_\e[92mi\e[39m\e[36ma\e[39m(δ_pq\e[36ma\e[39m g_\e[92mi\e[39m\e[36ma\e[39m\e[92mi\e[39m\e[92mi\e[39m h_p\e[36ma\e[39m E_pq)"
    @test string(SpinAdaptedSecondQuantization.get_all_indices(t)) == "SpinAdaptedSecondQuantization.MOIndex[p, q, \e[92mi\e[39m, \e[36ma\e[39m]"

    t = SpinAdaptedSecondQuantization.lower_delta_indices(t)

    @test string(t) == "3/5 ∑_\e[92mi\e[39m\e[36ma\e[39m(δ_pq\e[36ma\e[39m g_\e[92mi\e[39mp\e[92mi\e[39m\e[92mi\e[39m h_pp E_pp)"
end

@testset "expression constructors" begin
    p = general(1)
    q = general(2)
    i = occupied(1)
    a = virtual(1)

    Epq = E(p, q)
    dpi = δ(p, i)
    dai = δ(a, i)
    hai = real_tensor("h", a, i)

    @test string(Epq) == "E_pq"
    @test string(dpi) == "δ_p\e[92mi\e[39m"
    @test string(dai) == "0"
    @test string(hai) == "h_\e[36ma\e[39m\e[92mi\e[39m"
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

    t3 = SpinAdaptedSecondQuantization.make_space_for_index(t, i)
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

@testset "term summation delta" begin
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

    t = SpinAdaptedSecondQuantization.lower_delta_indices(t)

    #p -> a
    t1 = ∑(t, [p])
    @test t1 == ∑(δ(a,p) * real_tensor("h", p, i) * E(p, q), [p]).terms[1]
    t1 = SpinAdaptedSecondQuantization.simplify_summation_deltas(t1)
    @test t1 == (real_tensor("h", a, i) * E(a, q)).terms[1]
    t2 = ∑(t1, [a])
    @test t2 == ∑(real_tensor("h", a, i) * E(a, q), [a]).terms[1]

    #a -> p
    t1 = ∑(t, [a])
    @test t1 == ∑(δ(a,p) * real_tensor("h", p, i) * E(p, q), [a]).terms[1]
    t1 = SpinAdaptedSecondQuantization.simplify_summation_deltas(t1)
    @show t1
    t2 = ∑(t1, [p])
    @test t2 == ∑(real_tensor("h", a, i) * E(a, q), [a]).terms[1]

    #ap
    t1 = ∑(t, [a, p])
    @test t1 == ∑(δ(a,p) * real_tensor("h", p, i) * E(p, q), [p, a]).terms[1]
    t1 = SpinAdaptedSecondQuantization.simplify_summation_deltas(t1)
    @test t1 == ∑(real_tensor("h", a, i) * E(a, q), [a]).terms[1]

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

@testset "term multiplication" begin
    println()

    p = general(1)
    q = general(2)
    r = general(3)
    s = general(4)
    i = occupied(1)
    j = occupied(2)
    a = virtual(1)
    b = virtual(2)

    t1 = SpinAdaptedSecondQuantization.Term(
        3 // 5,
        SpinAdaptedSecondQuantization.MOIndex[i, j],
        SpinAdaptedSecondQuantization.KroneckerDelta[],
        SpinAdaptedSecondQuantization.Tensor[],
        [SpinAdaptedSecondQuantization.SingletExcitationOperator(i, j)],
    )

    @show t1

    t2 = t1 * t1

    @show t2

    println()

    t1 = SpinAdaptedSecondQuantization.Term(
        3 // 5,
        SpinAdaptedSecondQuantization.MOIndex[],
        [
            SpinAdaptedSecondQuantization.KroneckerDelta(i, p)
        ],
        SpinAdaptedSecondQuantization.Tensor[],
        [SpinAdaptedSecondQuantization.SingletExcitationOperator(p, q)],
    )

    t2 = SpinAdaptedSecondQuantization.Term(
        3 // 5,
        SpinAdaptedSecondQuantization.MOIndex[],
        [
            SpinAdaptedSecondQuantization.KroneckerDelta(a, p)
        ],
        SpinAdaptedSecondQuantization.Tensor[
            SpinAdaptedSecondQuantization.RealTensor("h", [p, q])
        ],
        SpinAdaptedSecondQuantization.Operator[],
    )

    @show t1 t2 t1 * t2

    println()
end

@testset "expression addition" begin
    p = general(1)
    q = general(2)
    r = general(3)
    s = general(4)
    i = occupied(1)
    j = occupied(2)
    a = virtual(1)
    b = virtual(2)

    e1 = δ(a, b)
    e2 = real_tensor("h", i, j)

    @test (e1 + 1) == (1 + e1)
    @test (e1 + e2) == (e2 + e1)
    @test (e1 + e2 + 2) == (2 + e1 + e2) == (e2 + 2 + e1)
end

@testset "expression multiplication" begin
    p = general(1)
    q = general(2)
    r = general(3)
    s = general(4)
    i = occupied(1)
    j = occupied(2)
    k = occupied(3)
    l = occupied(4)
    m = occupied(5)
    n = occupied(6)
    a = virtual(1)
    b = virtual(2)

    e1 = δ(a, b)
    e2 = real_tensor("h", i, j)
    h = ∑(real_tensor("h", i, j) * E(i, j), [i, j])
    g = ∑(real_tensor("g", i, j, k, l) * E(i, j) * E(k, l), [i, j, k, l])

    @test ((e1 + e2) * 3) // 5 * e2 == 3//5 * δ(a,b) * real_tensor("h", i, j) + 3//5 * real_tensor("h", i, j) * real_tensor("h", i, j)
    @test (h + h) == 2 * h
    @test simplify(g * h) == ∑(real_tensor("g", i, j, k, l) * real_tensor("h", m, n) * E(i, j) * E(k, l) * E(m, n), [i,j,k,l,m,n])
    @test simplify(h * g) == ∑(real_tensor("g", i, j, k, l) * real_tensor("h", m, n) * E(m, n) * E(i, j) * E(k, l), [i,j,k,l,m,n])
end

@testset "simple commutator" begin
    zero = SpinAdaptedSecondQuantization.Expression([SpinAdaptedSecondQuantization.Term(0,
        SpinAdaptedSecondQuantization.MOIndex[],
        SpinAdaptedSecondQuantization.KroneckerDelta[],
        SpinAdaptedSecondQuantization.Tensor[],
        SpinAdaptedSecondQuantization.SingletExcitationOperator[],
    )])

    p = general(1)
    q = general(2)
    r = general(3)
    s = general(4)
    i = occupied(1)
    j = occupied(2)
    k = occupied(3)
    l = occupied(4)
    a = virtual(1)
    b = virtual(2)

    @test commutator(E(p, q), E(r, s)) == E(p, s) * δ(q, r) - E(r, q) * δ(p, s)
    @test commutator(E(i, j), E(a, k)) == - E(a, j) * δ(i, k)
    @test commutator(E(i, a), E(a, i)) == E(i, i) - E(a, a)
    @test commutator(E(a, i), E(a, i)) == zero
    
    h = ∑(real_tensor("h", p, q) * E(p, q), [p, q])
    hcomm = simplify(commutator(h, E(p, q)))
    @test hcomm == summation(- real_tensor("h", q, r) * E(p, r) + real_tensor("h", r, p) * E(r, q), [r])
    @test hcomm * E(r, r) == summation(- real_tensor("h", q, s) * E(p, s) * E(r, r) + real_tensor("h", s, p) * E(s, q) * E(r, r), [s])
end
