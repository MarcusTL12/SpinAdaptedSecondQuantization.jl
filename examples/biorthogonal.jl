using SpinAdaptedSecondQuantization

h = ∑((
          real_tensor("F", 1, 2) +
          ∑((-2 // 1 * psym_tensor("g", 1, 2, 3, 3) +
             psym_tensor("g", 1, 3, 3, 2)) * occupied(3), [3])
      ) * E(1, 2) * electron(1, 2), 1:2)

g = 1 // 2 * ∑(psym_tensor("g", 1, 2, 3, 4) * e(1, 2, 3, 4) *
               electron(1, 2, 3, 4), 1:4)

H = simplify(h + g)

T2 = 1 // 2 * ∑(psym_tensor("t", 1, 2, 3, 4) * E(1, 2) * E(3, 4) *
                occupied(2, 4) * virtual(1, 3), 1:4)

# Compute full e^-T H e^T to truncation
Hbar = simplify(bch(H, T2, 4))

# Compute Hbar |HF⟩, keeping only terms with up to two operators
Hbar_ket = simplify(act_on_ket(Hbar, 2))

omega_ai = begin
    # Project ⟨ai| Hbar |HF⟩
    tmp = project_biorthogonal(Hbar_ket, E(1, 2) * occupied(2) * virtual(1))

    # Look for pattern: u_aibj = 2 t_aibj - t_ajbi
    look_for_tensor_replacements(tmp, make_exchange_transformer("t", "u"))
end

println("omega_ai = ", omega_ai)
println()

omega_aibj = begin
    # Project ⟨aibj| Hbar |HF⟩
    tmp = project_biorthogonal(Hbar_ket,
        E(1, 2) * E(3, 4) * occupied(2, 4) * virtual(1, 3))

    # Expand permutations (1, 2) <=> (3, 4)
    # to allow for better simplification
    tmp = symmetrize(tmp, make_permutation_mappings([(1, 2), (3, 4)]))
    tmp = simplify_heavy(tmp)

    # Look for pattern: u_aibj = 2 t_aibj - t_ajbi
    tmp = look_for_tensor_replacements(tmp, make_exchange_transformer("t", "u"))

    # Look for pairs of terms differing only by permuting (1, 2) <=> (3, 4)
    s, ss, _ = desymmetrize(tmp, make_permutation_mappings([(1, 2), (3, 4)]))

    # Combining self symmetric and symmetrizing terms
    # to one combined expression
    s + ss // 2
end

println("omega_aibj = P_ai,bj (", omega_aibj)
println(")\n")
