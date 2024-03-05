include("multilevel.jl")

right_state = fermiondag(1, α) * ivir(1)

function act_on_actual_ket(O)
    act_on_ket(O * right_state)
end

h_elec = ∑(real_tensor("h", 1, 2) * E(1, 2) * active(1, 2), 1:2)
h_posi = ∑(real_tensor("h", 1, 2) * E(1, 2) * ivir(1, 2), 1:2)

h = h_elec + h_posi

g_elec = 1 // 2 * ∑(psym_tensor("g", 1, 2, 3, 4) * e(1, 2, 3, 4) * active(1, 2, 3, 4), 1:4) |> simplify
g_posi = 1 // 2 * ∑(psym_tensor("g", 1, 2, 3, 4) * e(1, 2, 3, 4) * ivir(1, 2, 3, 4), 1:4) |> simplify

g_int = -∑(psym_tensor("g", 1, 2, 3, 4) * e(1, 2, 3, 4) * active(1, 2) * ivir(3, 4), 1:4) |> simplify

H_elec = h_elec + g_elec

H = H_elec + h_posi + g_posi + g_int

E_hf = hf_expectation_value(right_state' * H * right_state) |> simplify_heavy

@show E_hf

println()

hF_elec = ∑((real_tensor("F", 1, 2) +
             ∑(aocc(3) * (-2psym_tensor("g", 1, 2, 3, 3) +
                          psym_tensor("g", 1, 3, 3, 2)), [3])) * E(1, 2) * active(1, 2), 1:2)

HF_elec = simplify(hF_elec + g_elec)

HF = HF_elec + h_posi + g_posi + g_int

E_hf = hf_expectation_value(right_state' * (H + HF) // 2 * right_state) |> simplify_heavy

@show E_hf

println()

ex_ketop(a, i) = E(a, i) * aocc(i) * avir(a)
ex_ketop(a, i, b, j) = E(a, i) * E(b, j) * aocc(i, j) * avir(a, b)

deex_braop(a, i) = 1 // 2 * ex_ketop(a, i)'
deex_braop(a, i, b, j) = 1 // 3 * ex_ketop(a, i, b, j)' +
                         1 // 6 * ex_ketop(a, j, b, i)'

ex_positron(a, b) = E(a, b) * ivir(a, b)
ex_positron(a, b, c, d) = E(a, i) * E(b, j) * ivir(a, b, c, d)

deex_positron(a, i) = 1 // 2 * ex_positron(a, i)'
deex_positron(a, i, b, j) = 1 // 3 * ex_positron(a, i, b, j)' +
                            1 // 6 * ex_positron(a, j, b, i)'

t(inds...) = psym_tensor("t", inds...)
s(inds...) = psym_tensor("s", inds...)

T2 = 1 // 2 * ∑(
    t(1:4...) * ex_ketop(1, 2, 3, 4),
    1:4
)

S1 = ∑(s(1, 2, 3, 4) * ex_ketop(1, 2) * ex_positron(3, 4), 1:4)
S2 = 1 // 2 * ∑(
    s(1:6...) * ex_ketop(1, 2, 3, 4) * ex_positron(5, 6),
    1:6
)

T = T2 + S1 + S2

function omega(proj, op, n)
    hf_expectation_value(simplify(right_state' * proj * bch(op, T, n) * right_state))
end

function omega_AB()
    o = omega(deex_positron(2, 1), HF_elec, 1)
    o = simplify_heavy(o)
    o = look_for_tensor_replacements_smart(o, make_exchange_transformer("t", "u"))
    o = look_for_tensor_replacements_smart(o, make_exchange_transformer("g", "L"))
end

;
