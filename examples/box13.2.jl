using SpinAdaptedSecondQuantization

h = ∑((
        real_tensor("F", 1, 2) +
        ∑((-2 // 1 * psym_tensor("g", 1, 2, 3, 3) +
           psym_tensor("g", 1, 3, 3, 2)) * occupied(3), [3])
    ) * E(1, 2) * electron(1, 2), 1:2)

g = 1 // 2 * ∑(psym_tensor("g", 1, 2, 3, 4) * e(1, 2, 3, 4) * electron(1, 2, 3, 4), 1:4)

H = simplify(h + g)

Eai(a, i) = E(a, i) * virtual(a) * occupied(i)

# H |HF⟩
function eq0()
    simplify(act_on_ket(H))
end

# [H, Eai] |HF⟩
function eq1()
    x = simplify(act_on_ket(commutator(H, Eai(1, 2))))

    x = look_for_tensor_replacements_smart(x, make_exchange_transformer("g", "L"))

    x
end

# [[H, Eai], Ebj] |HF⟩
function eq2()
    x = simplify(act_on_ket(
        commutator(commutator(H, Eai(1, 2)), Eai(3, 4))
    ))

    x = look_for_tensor_replacements_smart(x, make_exchange_transformer("g", "L"))

    s, ss, ns = desymmetrize(x, make_permutation_mappings([(1, 2), (3, 4)]))

    println("Symmetrize:")
    println(s)

    println("\nSymmetric:")
    println(ss)

    println("\nNon-symmetric:")
    println(ns)
end

# [[[H, Eai], Ebj], Eck] |HF⟩
function eq3()
    x = simplify(act_on_ket(
        commutator(commutator(commutator(H, Eai(1, 2)), Eai(3, 4)), Eai(5, 6))
    ))

    x = look_for_tensor_replacements_smart(x, make_exchange_transformer("g", "L"))

    s, ss, ns = desymmetrize(x, make_permutation_mappings([(1, 2), (3, 4), (5, 6)]))

    println("Symmetrize:")
    println(s)

    println("\nSymmetric:")
    println(ss)

    println("\nNon-symmetric:")
    println(ns)
end

# [[[[H, Eai], Ebj], Eck], Edl] |HF⟩
function eq4()
    x = simplify(act_on_ket(
        commutator(commutator(commutator(commutator(H, Eai(1, 2)), Eai(3, 4)), Eai(5, 6)), Eai(7, 8))
    ))

    x = look_for_tensor_replacements_smart(x, make_exchange_transformer("g", "L"))

    s, ss, ns = desymmetrize(x, make_permutation_mappings([(1, 2), (3, 4), (5, 6), (7, 8)]))

    println("Symmetrize:")
    println(s)

    println("\nSymmetric:")
    println(ss)

    println("\nNon-symmetric:")
    println(ns)
end
