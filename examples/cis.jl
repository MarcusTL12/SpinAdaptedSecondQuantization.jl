using SpinAdaptedSecondQuantization

hF = ∑(
    (rsym_tensor("F", 1, 2) + ∑(
        (-2rsym_tensor("g", 1, 2, 3, 3) + rsym_tensor("g", 1, 3, 3, 2)) *
        occupied(3),
        [3]
    )) * E(1, 2),
    1:2
)
g = 1 // 2 * ∑(rsym_tensor("g", 1:4...) * e(1:4...), 1:4) |> simplify

HF = simplify(hF + g)

EHF = act_on_ket(HF, 0) |> simplify

E_single = simplify(
    (act_on_ket(E(1, 2) * HF * E(2, 1), 0) // 2 -
     act_on_ket(HF, 0)) * occupied(1) * virtual(2)
)

@show E_single
println()

H_cis = simplify_heavy(
    (act_on_ket(E(2, 3) * HF * E(4, 1), 0) // 2 -
     act_on_ket(HF, 0) * δ(1, 2) * δ(3, 4)) *
    occupied(1, 2) * virtual(3, 4)
)

@show H_cis
