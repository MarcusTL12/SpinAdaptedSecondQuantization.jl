using Test
using SpinAdaptedSecondQuantization
import SpinAdaptedSecondQuantization as SASQ

SASQ.enable_color()
@testset "ccsd equations" begin
    Φ = 1 // 2 * ∑(psym_tensor("g", 1,2,3,4) * e(1,2,3,4), 1:4) +
        ∑((-2 * psym_tensor("g", 1, 2, 3, 3) + psym_tensor("g", 1, 3, 3, 2)) *
            occupied(3) * E(1,2), 1:3)
    F = ∑(real_tensor("F", 1, 2) * E(1, 2), 1:2)
    E_ref = ∑(2 * real_tensor("F", 1, 1) * occupied(1), [1]) +
            ∑(-2* psym_tensor("g", 1, 1, 2, 2) * occupied(1,2), [1,2]) +
            ∑(1 * psym_tensor("g", 1, 2, 2, 1) * occupied(1,2), [1,2])

    H = F + Φ - E_ref |> simplify

    Tn(n) = 1 // factorial(n) * ∑(psym_tensor("t", 1:2n...) * prod(E(2i-1,2i) for i = 1:n) * occupied(2:2:2n...) * virtual(1:2:2n...), 1:2n)
    T2 = Tn(2)
    T3 = Tn(3)
    T4 = Tn(4)

    P1 = 1//2 * E(1,2) * occupied(1) * virtual(2)
    P2 = 1//6 * (2*E(1,2) * E(3,4) + E(3,2) * E(1,4)) * occupied(1,3) * virtual(2,4)


    HT2 = bch(H, T2, 2) |> act_on_ket |> simplify

    E_terms = [length(t.operators) == 0 for t in HT2.terms]
    Ecorr = SASQ.Expression(HT2[E_terms]) |> simplify_heavy
    @show Ecorr

    ai_terms = [length(t.operators) == 1 for t in HT2.terms]
    omega_ai = SASQ.Expression(HT2[ai_terms]) |> simplify
    omega_ai = P1 * omega_ai |> act_on_ket |> simplify_heavy
    @show omega_ai

    aibj_terms = [length(t.operators) == 2 for t in HT2.terms]
    omega_aibj = SASQ.Expression(HT2[aibj_terms]) |> simplify
    omega_aibj = P2 * omega_aibj |> act_on_ket |> simplify_heavy
    @show omega_aibj
end
SASQ.disable_color()