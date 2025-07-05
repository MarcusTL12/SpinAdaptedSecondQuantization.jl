include("multilevel.jl")

h = ∑(real_tensor("h", 1, 2) * E(1, 2) * electron(1, 2), 1:2)
g = 1 // 2 * ∑(psym_tensor("g", 1:4...) * e(1:4...) * electron(1:4...), 1:4) |> simplify

H = simplify(h + g)

hF = ∑((real_tensor("F", 1, 2) +
        ∑(occupied(3) * (-2psym_tensor("g", 1, 2, 3, 3) +
                         psym_tensor("g", 1, 3, 3, 2)), [3])) * E(1, 2) * electron(1, 2), 1:2)

HF = simplify(hF + g)

EHF = hf_expectation_value((H + HF) // 2) |> simplify_heavy

ex_ketop(a, i) = E(a, i) * aocc(i) * avir(a)
ex_ketop(a, i, b, j) = E(a, i) * E(b, j) * aocc(i, j) * avir(a, b)

deex_braop(a, i) = 1 // 2 * ex_ketop(a, i)'
deex_braop(a, i, b, j) = 1 // 3 * ex_ketop(a, i, b, j)' +
                         1 // 6 * ex_ketop(a, j, b, i)'

t(inds...) = psym_tensor("t", inds...)

T1 = ∑(t(1, 2) * ex_ketop(1, 2), 1:2)
T2 = 1 // 2 * ∑(
    t(1:4...) * ex_ketop(1, 2, 3, 4),
    1:4
)

function get_E_corr_t1()
    E = hf_expectation_value(commutator(HF, T2)) |> simplify_heavy
    E = look_for_tensor_replacements_smart(E, make_exchange_transformer("t", "u"))
    E = look_for_tensor_replacements_smart(E, make_exchange_transformer("g", "L"))
end

function omega(proj, op, n)
    hf_expectation_value(simplify(proj * bch(op, T2, n)))
end

function get_omega_t1()
    o = omega(deex_braop(1, 2), HF, 2)
    o = simplify_heavy(o)
    o = look_for_tensor_replacements_smart(o, make_exchange_transformer("t", "u"))
    o = look_for_tensor_replacements_smart(o, make_exchange_transformer("g", "L"))
end

function get_omega_t2()
    o = omega(deex_braop(1, 2, 3, 4), HF, 2)
    o = simplify_heavy(o)
    o = look_for_tensor_replacements_smart(o, make_exchange_transformer("t", "u"))
    o = look_for_tensor_replacements_smart(o, make_exchange_transformer("g", "L"))

    s, ss, _ = desymmetrize(o, make_permutation_mappings([(1, 2), (3, 4)]))

    println("s: $s\n")
    println("ss: $ss")
end
