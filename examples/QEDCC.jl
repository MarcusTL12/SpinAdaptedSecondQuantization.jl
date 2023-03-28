using SpinAdaptedSecondQuantization

function cc_projection(order)
    # Projection operators
    El(i,a) = E(i,a) * virtual(a)  * occupied(i)
    if order == 0
        P = SASQ.Expression(1)
    elseif order == 1
        P = 1//2 * El(1,2)
    elseif order == 2
        P = 1//6 * (2 * El(1,2) * El(3,4) + El(3,2) * El(1,4))
    else
        throw("order not supported")
    end
    return P
end

function cc_ket(H, T, n, order)
    # Choose projection
    P = cc_projection(order)
    # HT = e^-T H eT |HF>
    HT = act_on_ket(bch(H, T, n), order) |> simplify

    # Return only terms of op_length = order
    terms = [length(t.operators) == order for t in HT.terms]
    ex = P * SASQ.Expression(HT[terms]) |> act_on_ket |> simplify_heavy
    return ex
end

function project(HT, order)
    terms = HT[[length(t.operators) == order for t in HT.terms]]
    if order > 1
        for (i, t) = enumerate(terms)
            terms[i] = SASQ.Term(
            t.scalar,
            t.sum_indices,
            t.deltas,
            vcat(t.tensors, SASQ.ParticleSymmetricTensor("A", collect(1:2*order))),
            SASQ.Operator[],
            t.constraints
            )
        end
        ex = SASQ.Expression(terms)
        terms = simplify_heavy(ex).terms
    end
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
    return simplify_heavy(SASQ.Expression(terms))
end

function cc_ket2(H, T, n, order)
    # HT = e^-T H eT |HF>
    HT = bch(H, T, n) |> act_on_ket |> simplify

    # Return only terms of op_length = order
    return project(HT, order)
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
Γ = real_tensor("γ") * bosondag()

# HT = e^-T H eT |HF>
H = F + Φ + d + E_f |> simplify
T = Tn(2) + Sn(1) + Sn(2) + Γ # + Tn(3) + Tn(4)
HT = act_on_ket(bch(H, T, 2)) |> simplify;

order = 0
P = cc_projection(order)
terms = [sum(typeof.(t.operators) .== SASQ.SingletExcitationOperator) == order for t in HT.terms];
ex = P * SASQ.Expression(HT[terms]);
ex = ex |> act_on_ket |> simplify_heavy;
phorder = 0
terms = [sum(typeof.(t.operators) .== SASQ.BosonOperator) == phorder for t in ex.terms];
Ecorr = SASQ.Expression(ex[terms])

phorder = 1
terms = [sum(typeof.(t.operators) .== SASQ.BosonOperator) == phorder for t in ex.terms];
omega_1 = SASQ.Expression(ex[terms])

order = 1
P = cc_projection(order)
terms = [sum(typeof.(t.operators) .== SASQ.SingletExcitationOperator) == order for t in HT.terms];
ex = P * simplify_heavy(SASQ.Expression(HT[terms]));
ex = ex |> act_on_ket |> simplify_heavy;
phorder = 0
terms = [sum(typeof.(t.operators) .== SASQ.BosonOperator) == phorder for t in ex.terms];
omega_ai = SASQ.Expression(ex[terms])

phorder = 1
terms = [sum(typeof.(t.operators) .== SASQ.BosonOperator) == phorder for t in ex.terms];
omega_ai_1 = SASQ.Expression(ex[terms])

order = 2
P = cc_projection(order)
terms = [sum(typeof.(t.operators) .== SASQ.SingletExcitationOperator) == order for t in HT.terms];
ex = P * simplify_heavy(SASQ.Expression(HT[terms]))
ex = ex |> act_on_ket |> simplify_heavy;
phorder = 0
terms = [sum(typeof.(t.operators) .== SASQ.BosonOperator) == phorder for t in ex.terms];
omega_aibj = SASQ.Expression(ex[terms])

phorder = 1
terms = [sum(typeof.(t.operators) .== SASQ.BosonOperator) == phorder for t in ex.terms];
omega_aibj_1 = SASQ.Expression(ex[terms]) |> simplify_heavy