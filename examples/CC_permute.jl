using SpinAdaptedSecondQuantization

function project(HT, order)
    # <aibj|ckdl> = P^ab_ij δ_ab δ_ij
    # Replace all operators E_ai E_bj E_ck with P_aibjck and remove summation
    # Assumes all operators are of form E_ai and sorted and summed over
    terms = HT[[length(t.operators) == order for t in HT.terms]]
    if order > 1
        for (i, t) = enumerate(terms)
            terms[i] = SASQ.Term(
            t.scalar,
            t.sum_indices[2*order+1:end],
            t.deltas,
            vcat(t.tensors, SASQ.PermuteTensor(collect(1:2*order))),
            SASQ.Operator[],
            t.constraints
            )
        end
    else
        for (i, t) = enumerate(terms)
            terms[i] = SASQ.Term(
            t.scalar,
            t.sum_indices[2*order+1:end],
            t.deltas,
            t.tensors,
            SASQ.Operator[],
            t.constraints
            )
        end
    end
    return simplify_heavy(SASQ.Expression(terms))
end

function cc_ket_P(H, T, n, order)
    # HT = e^-T H eT |HF>
    HT = bch(H, T, n) |> x -> act_on_ket(x, order) |> simplify

    # Return only terms of op_length = order
    return project(HT, order)
end

# Electronic cluster operator to n'th order
# T = 1/N! ∑_μ t_μ τ_μ
τ(n) = prod(E(2i-1,2i) for i = 1:n) * occupied(2:2:2n...) * virtual(1:2:2n...)
Tn(n) = 1 // factorial(n) * ∑(psym_tensor("t", 1:2n...) * τ(n), 1:2n)

# Define Hamiltonian in terms of F and g
Φ = 1 // 2 * ∑(psym_tensor("g", 1,2,3,4) * e(1,2,3,4), 1:4) +
    ∑((-2 * psym_tensor("g", 1, 2, 3, 3) + psym_tensor("g", 1, 3, 3, 2)) *
        occupied(3) * E(1,2), 1:3)
F = ∑(real_tensor("F", 1, 2) * E(1, 2), 1:2)
H = F + Φ |> simplify

# Solve CCSDTQ equations (~2 min / 8 threads)
# Note: P_aibj is the Permutation operator
T = Tn(2) + Tn(3) + Tn(4)
Ecorr = cc_ket_P(H, T, 1, 0)
omega_ai = cc_ket_P(H, T, 2, 1)
@time omega_aibj = cc_ket_P(H, T, 2, 2)
@time omega_aibjck = cc_ket_P(H, T, 2, 3)
@time omega_aibjckdl = cc_ket_P(H, T, 3, 4)