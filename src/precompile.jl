using PrecompileTools

@compile_workload begin
    # Covering Epq and realtensor
    h = ∑(real_tensor("h", 1, 2) * E(1, 2), 1:2)
    g = 1 // 2 * ∑(real_tensor("g", 1:4...) * e(1:4...), 1:4) |> simplify

    H = simplify(h + g)

    hF = ∑(
        (real_tensor("F", 1, 2) + ∑(
            (-2real_tensor("g", 1, 2, 3, 3) + real_tensor("g", 1, 3, 3, 2)) *
            occupied(3),
            [3]
        )) * E(1, 2),
        1:2
    )
    HF = hF + g

    simplify_heavy(act_on_ket((HF + H) // 2, 0))

    # rsym_tensor and psym_tensor
    h = ∑(rsym_tensor("h", 1, 2) * E(1, 2), 1:2)
    g = 1 // 2 * ∑(psym_tensor("g", 1:4...) * e(1:4...), 1:4) |> simplify

    H = simplify(h + g)

    hF = ∑(
        (rsym_tensor("F", 1, 2) + ∑(
            (-psym_tensor("g", 1, 2, 3, 3) + psym_tensor("g", 1, 3, 3, 2)) *
            occupied(3),
            [3]
        )) * E(1, 2),
        1:2
    )
    HF = hF + g

    simplify_heavy(act_on_ket((HF + H) // 2, 0))

    # Add more simple calls to precompile other types/functions
end
