using SpinAdaptedSecondQuantization
using Permutations


function project(HT, order, photon_order)
    # <aibj|ckdl> = P^ab_ij δ_ab δ_ij
    # Replace all operators E_ai E_bj E_ck with P_aibjck and remove summation
    # Assumes all operators are of form E_ai and sorted and summed over
    n_E = [sum([typeof(op)==SASQ.SingletExcitationOperator for op in t.operators], init=0) == order for t in HT.terms]
    n_bosons = [sum([typeof(op)==SASQ.BosonOperator for op in t.operators], init=0) == photon_order for t in HT.terms]
    terms = HT[n_E .* n_bosons]

    for (i, t) = enumerate(terms)
        if order > 1
            new_tensors = vcat(t.tensors, SASQ.PermuteTensor(collect(1:2*order)))
        else
            new_tensors = t.tensors
        end
        terms[i] = SASQ.Term(
        t.scalar,
        t.sum_indices[2*order+1:end],
        t.deltas,
        new_tensors,
        SASQ.Operator[],
        t.constraints
        )
    end
    return simplify_heavy(SASQ.Expression(terms))
end


function cc_ket_P(H, T, n, order, photon_order)
    # HT = e^-T H eT |HF>
    HT = bch(H, T, n) |> x -> act_on_ket(x, order + photon_order) |> simplify

    # Return only terms of op_length = order
    return project(HT, order, photon_order)
end

function remove_P(expression)
    ex = SASQ.Expression(deepcopy(expression.terms))
    for term in copy(ex.terms)
        for (i, tensor) in enumerate(term.tensors)
            if typeof(tensor) == SASQ.PermuteTensor
                deleteat!(term.tensors, i)
                break
            end
        end
    end
    return ex
end


function expand_P(expression)
    ex = SASQ.Expression(deepcopy(expression.terms))
    new_terms = deepcopy(ex.terms)
    empty!(new_terms)
    for term in copy(ex.terms)
        for (i, tensor) in enumerate(term.tensors)
            if typeof(tensor) == SASQ.PermuteTensor
                occ = tensor.indices[2:2:end]
                vir = tensor.indices[1:2:end]
                deleteat!(term.tensors, i)

                permutations = PermGen(length(occ))
                for p in permutations
                    new_term = SASQ.exchange_indices(deepcopy(term), vcat([vir[i] => vir[j] for (i,j) in enumerate(p.data)], 
                                                                          [occ[i] => occ[j] for (i,j) in enumerate(p.data)]))
                    push!(new_terms, new_term)
                end

                break
            end
        end
    end
    return SASQ.Expression(new_terms)
end


# Define Hamiltonian in terms of F and g
Φ = 1 // 2 * ∑(psym_tensor("g", 1,2,3,4) * e(1,2,3,4), 1:4) +
    ∑((-2 * psym_tensor("g", 1, 2, 3, 3) + psym_tensor("g", 1, 3, 3, 2)) *
        occupied(3) * E(1,2), 1:3)
F = ∑(real_tensor("F", 1, 2) * E(1, 2), 1:2)
d = ∑(real_tensor("d", 1,2) * (boson() + bosondag()) * E(1,2), 1:2)
E_f = real_tensor("ω") * bosondag() * boson()

# Electronic cluster operator to n'th order
Tn(n) = 1 // factorial(n) * ∑(psym_tensor("t", 1:2n...) * prod(E(2i-1,2i) for i = 1:n) * occupied(2:2:2n...) * virtual(1:2:2n...), 1:2n)
Sn(n) = 1 // factorial(n) * ∑(psym_tensor("s", 1:2n...) * bosondag() * prod(E(2i-1,2i) for i = 1:n) * occupied(2:2:2n...) * virtual(1:2:2n...), 1:2n)
T2u = ∑(1 // 2 * occupied(1, 3) * virtual(2, 4) *
    (2 // 3 * psym_tensor("u", 1, 2, 3, 4) +
    1 // 3 * psym_tensor("u", 1, 4, 3, 2)) * E(1, 2) * E(3, 4), 1:4)

# HT = e^-T H eT |HF>
gamma_transform = real_tensor("γ") * real_tensor("ω") * bosondag() +
                  ∑(real_tensor("d", 1,2) * real_tensor("γ") * E(1,2), 1:2)
H0 = F + Φ
H = H0 + d + E_f + gamma_transform
T = Tn(2) + Sn(1) + Sn(2) # Γ1 included directly in H

# New terms compared to electronic Hamiltonian CCSD
Ecorr = cc_ket_P(H, T, 1, 0, 0) - cc_ket_P(H0, Tn(2), 2, 0, 0)
omega_ai = cc_ket_P(H, T, 2, 1, 0) - cc_ket_P(H0, Tn(2), 2, 1, 0)
omega_aibj = cc_ket_P(H, T, 2, 2, 0) - cc_ket_P(H0, Tn(2), 2, 2, 0)
omega_1 = cc_ket_P(H, T, 1, 0, 1)
omega_ai_1 = cc_ket_P(H, T, 2, 1, 1)
omega_aibj_1 = cc_ket_P(H, T, 2, 2, 1)