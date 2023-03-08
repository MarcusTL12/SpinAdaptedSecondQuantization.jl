using Test
using SpinAdaptedSecondQuantization
import SpinAdaptedSecondQuantization as SASQ
SASQ.disable_color()

@testset "ccsd equations" begin
    Φ = 1 // 2 * ∑(psym_tensor("g", 1,2,3,4) * e(1,2,3,4), 1:4) +
        ∑((-2 * psym_tensor("g", 1, 2, 3, 3) + psym_tensor("g", 1, 3, 3, 2)) *
            occupied(3) * E(1,2), 1:3)
    F = ∑(real_tensor("F", 1, 2) * E(1, 2), 1:2)
    E_ref = ∑(2 * real_tensor("F", 1, 1) * occupied(1), [1]) +
            ∑(-2* psym_tensor("g", 1, 1, 2, 2) * occupied(1,2), [1,2]) +
            ∑(1 * psym_tensor("g", 1, 2, 2, 1) * occupied(1,2), [1,2])

    H = F + Φ - E_ref
    Hket = act_on_ket(H) |> simplify
    @show Hket

    T2 = 1 // 2 * ∑(psym_tensor("t", 1,2,3,4) * E(1,2)*E(3,4) * occupied(2,4)*virtual(1,3), 1:4)
    T3 = 1 // 6 * ∑(psym_tensor("t", 1,2,3,4,5,6) * E(1,2)*E(3,4)*E(5,6) * occupied(2,4,6)*virtual(1,3,5), 1:6)
    @show typeof(H), typeof(T2)
    @show T2
    @show T3
    @time HT2 = act_on_ket(commutator(H,T2)) |> simplify
    @time HT22 = 1//2 * act_on_ket(H * T2 * T2) |> simplify
    @time HT3 = act_on_ket(H * T3) |> simplify
    @time HT23 = act_on_ket(H * T3) |> simplify

    @show HT2
end