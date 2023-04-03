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

EHF = hf_expectation_value(HF) |> simplify

E_single = simplify(
    (hf_expectation_value(E(1, 2) * HF * E(2, 1)) // 2 -
     hf_expectation_value(HF)) * occupied(1) * virtual(2)
)

@show E_single

H_cis = simplify_heavy(
    (hf_expectation_value(E(2, 3) * HF * E(4, 1)) // 2 -
     hf_expectation_value(HF) * δ(1, 2) * δ(3, 4)) *
    occupied(1, 2) * virtual(3, 4)
)

@show H_cis
