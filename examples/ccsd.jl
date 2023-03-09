using SpinAdaptedSecondQuantization

h = ∑(rsym_tensor("h", 1, 2) * E(1, 2), 1:2)
g = 1 // 2 * ∑(rsym_tensor("g", 1:4...) * e(1:4...), 1:4) |> simplify
H = h + g

hF = ∑((rsym_tensor("F", 1, 2) +
        ∑(occupied(3) * (-2rsym_tensor("g", 1, 2, 3, 3) +
                         rsym_tensor("g", 1, 3, 3, 2)), [3])) * E(1, 2), 1:2)

HF = simplify(hF + g)

@show H HF

EHF = hf_expectation_value((H + HF) // 2) |> simplify

@show EHF

deex_braop(a, i, b, j) =
    (1 // 3 * E(i, a) * E(j, b) + 1 // 6 * E(i, b) * E(j, a))

T2 = ∑(1 // 2 * occupied(1, 3) * virtual(2, 4) *
       psym_tensor("t", 1, 2, 3, 4) * E(1, 2) * E(3, 4), 1:4)

T2u = ∑(1 // 2 * occupied(1, 3) * virtual(2, 4) *
        (2 // 3 * psym_tensor("u", 1, 2, 3, 4) +
         1 // 3 * psym_tensor("u", 1, 4, 3, 2)) * E(1, 2) * E(3, 4), 1:4)

@show T2 T2u

E_corr = hf_expectation_value(commutator(H, T2u)) |> simplify_heavy

@show E_corr

;
