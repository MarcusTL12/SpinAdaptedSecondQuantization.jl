# Timings for Deriving Coupled Cluster Equations

The benchmarking script below can be used as a test of performance
for deriving coupled cluster ground state equations for arbitrary order.

Here are the timings of this being run for a few truncation orders running with
48 threads on a dual `Intel(R) Xeon(R) Gold 6342 CPU @ 2.80GHz` cpu system

| Derivation Step | CCS     | CCSD   | CCSDT | CCSDTQ | CCSDTQP | CCSDTQP6 |
| --------------- | ------- | ------ | ----- | ------ | ------- | -------- |
| bch             | 17.1 µs | 58 ms  | 1.0 s |  8.3 s | 49.4 s  |  210 s   |
| simplify        | 0.2 ms  | 39 ms  | 0.6 s |  5.5 s | 28.7 s  |  127 s   |
| act_on_ket      | 0.3 ms  | 27 ms  | 1.3 s | 10.4 s | 74.2 s  |  423 s   |
| simplify        | 0.3 ms  | 28 ms  | 0.5 s |  5.0 s | 29.5 s  |  130 s   |
| finalize        | 0.3 ms  | 25 ms  | 0.1 s |  0.7 s |  7.7 s  |  113 s   |
| total           | 1.1 ms  | 0.2 s  | 3.5 s | 30.1 s |  190 s  | 1002 s   |

```julia
using SpinAdaptedSecondQuantization
using BenchmarkTools

h = ∑((
          real_tensor("F", 1, 2) +
          ∑((-2 * psym_tensor("g", 1, 2, 3, 3) +
             psym_tensor("g", 1, 3, 3, 2)) * occupied(3), [3])
      ) * E(1, 2) * electron(1, 2), 1:2)

g = 1 // 2 * simplify(
    ∑(psym_tensor("g", 1:4...) * e(1:4...) * electron(1:4...), 1:4)
)

H = h + g

Eai(a, i) = E(a, i) * virtual(a) * occupied(i)
Eai(a, i, rest...) = Eai(a, i) * Eai(rest...)

function make_T(n)
    1 // factorial(n) * ∑(psym_tensor("t", 1:2n...) * Eai(1:2n...), 1:2n)
end

exc_trans_t2u = make_exchange_transformer("t", "u")
exc_trans_g2L = make_exchange_transformer("g", "L")

function tensor_replacements(ex)
    ex = look_for_tensor_replacements(ex, exc_trans_t2u)
    look_for_tensor_replacements(ex, exc_trans_g2L)
end

function simplify_and_extract_omega_blocks(Hbar_ket_simplified, ::Val{N}) where N
    omegas = NTuple{3,SASQ.Expression{Rational{Int}}}[]

    for n in 1:N
        template_ket = Eai(1:2n...)

        omega_nonsym = project_biorthogonal(Hbar_ket_simplified, template_ket)

        if n == 1
            omega_final = (omega_nonsym,
                zero(SASQ.Expression{Rational{Int}}),
                zero(SASQ.Expression{Rational{Int}}))
        else
            perm_maps = make_permutation_mappings([(i, i + 1) for i in 1:2:2n])

            omega_sym = symmetrize(omega_nonsym, perm_maps)

            omega_simplified = simplify_heavy(omega_sym)

            omega_replaced = tensor_replacements(omega_simplified)

            omega_final = desymmetrize(omega_replaced, perm_maps)
        end

        push!(omegas, omega_final)
    end

    omegas
end

num_equals = 60

function run_benchmarks(valN::Val{N}) where {N}
    tot_min_time = 0.0
    tot_median_time = 0.0

    Ts = [make_T(n) for n in 2:N]

    println("\n\n\n" * "="^num_equals)
    println("New Benchmark for N = ", N, ":")
    println("Running with ", Threads.nthreads(), " threads")
    println("\n\n")

    println("\nMaking Hbar")
    b = @benchmark global Hbar = bch($H, $Ts, 4)
    display(b)
    tot_min_time += minimum(b).time
    tot_median_time += median(b).time

    println("\nSimplifying Hbar")
    b = @benchmark global Hbar_simplified = simplify($Hbar)
    display(b)
    tot_min_time += minimum(b).time
    tot_median_time += median(b).time

    println("\nActing Hbar on |HF⟩")
    b = @benchmark global Hbar_ket = act_on_ket($Hbar_simplified)
    display(b)
    tot_min_time += minimum(b).time
    tot_median_time += median(b).time

    println("\nSimplifying Hbar_ket")
    b = @benchmark global Hbar_ket_simplified = simplify($Hbar_ket)
    display(b)
    tot_min_time += minimum(b).time
    tot_median_time += median(b).time

    println("\nBiorthonormal projections and simplifications")
    b = @benchmark global omegas = simplify_and_extract_omega_blocks(
        $Hbar_ket_simplified, $valN)
    display(b)
    tot_min_time += minimum(b).time
    tot_median_time += median(b).time

    print("\n\nAccumulate minimum benchmark times = ",
        tot_min_time * 1e-9)
    println(" s\n")

    print("Accumulate median benchmark times = ",
        tot_median_time * 1e-9)
    println(" s\n")

    nterms_tot = 0

    println("Number of terms listed by order:")

    for (n, (r, s, ns)) in enumerate(omegas)
        nterms = 0
        iszero(r) || (nterms += length(r.terms))
        iszero(s) || (nterms += length(s.terms))
        iszero(ns) || (nterms += length(ns.terms))

        nterms_tot += nterms

        println("n = $n => ", nterms, " terms")
    end

    println("\nTotal number of terms: ", nterms_tot)

    println("\n" * "="^num_equals)
end
```

Timing CCSDTQ equations as an example:

```julia
julia> run_benchmarks(Val(4))
============================================================
New Benchmark for N = 4:
Running with 24 threads




Making Hbar
BenchmarkTools.Trial: 1 sample with 1 evaluation per sample.
 Single result which took 7.030 s (24.51% GC) to evaluate,
 with a memory estimate of 26.45 GiB, over 405713131 allocations.

Simplifying Hbar
BenchmarkTools.Trial: 2 samples with 1 evaluation per sample.
 Range (min … max):  4.858 s …   4.992 s  ┊ GC (min … max): 11.23% … 15.90%
 Time  (median):     4.925 s              ┊ GC (median):    13.60%
 Time  (mean ± σ):   4.925 s ± 95.323 ms  ┊ GC (mean ± σ):  13.60% ±  3.30%

  █                                                       █  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  4.86 s         Histogram: frequency by time        4.99 s <

 Memory estimate: 6.27 GiB, allocs estimate: 100915447.

Acting Hbar on |HF⟩
BenchmarkTools.Trial: 1 sample with 1 evaluation per sample.
 Single result which took 13.359 s (44.81% GC) to evaluate,
 with a memory estimate of 100.77 GiB, over 1516750248 allocations.

Simplifying Hbar_ket
BenchmarkTools.Trial: 1 sample with 1 evaluation per sample.
 Single result which took 5.549 s (9.04% GC) to evaluate,
 with a memory estimate of 4.99 GiB, over 80958614 allocations.

Biorthonormal projections and simplifications
BenchmarkTools.Trial: 2 samples with 1 evaluation per sample.
 Range (min … max):  4.523 s …   4.625 s  ┊ GC (min … max): 37.86% … 36.09%
 Time  (median):     4.574 s              ┊ GC (median):    36.96%
 Time  (mean ± σ):   4.574 s ± 71.648 ms  ┊ GC (mean ± σ):  36.96% ±  1.25%

  █                                                       █  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  4.52 s         Histogram: frequency by time        4.62 s <

 Memory estimate: 11.95 GiB, allocs estimate: 129928374.


Accumulate minimum benchmark times = 35.318716928 s

Accumulate median benchmark times = 35.4367837025 s

Number of terms listed by order:
n = 1 => 11 terms
n = 2 => 24 terms
n = 3 => 36 terms
n = 4 => 57 terms

Total number of terms: 128

============================================================
```
