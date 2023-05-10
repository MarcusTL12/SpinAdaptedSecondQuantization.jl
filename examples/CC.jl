using SpinAdaptedSecondQuantization

function cc_projection(order)
    # Biorthogonal projection operators
    El(i,a) = E(i,a) * virtual(a) * occupied(i)
    a = 1; b = 3; c = 5
    i = 2; j = 4; k = 6
    if order == 0
        P = SASQ.Expression(1)
    elseif order == 1
        P = 1//2 * El(i,a)
    elseif order == 2
        P = 1//6 * (2 * El(i,a) * El(j,b) + El(i,b) * El(j,a))
    elseif order == 3
        P = 1//120 * (17 * El(i, a) * El(j, b) * El(k, c)
                     - 1 * El(i, a) * El(j, c) * El(k, b)
                     - 1 * El(i, b) * El(j, a) * El(k, c)
                     - 7 * El(i, b) * El(j, c) * El(k, a)
                     - 7 * El(i, c) * El(j, a) * El(k, b)
                     - 1 * El(i, c) * El(j, b) * El(k, a) )
    else
        throw("order not supported")
    end
    return P
end

function cc_ket(H, T, n, order)
    # <order| bch(H, T, n) | HF >

    P = cc_projection(order)
    HT = bch(H, T, n) |> act_on_ket |> simplify
    # Return only terms of op_length = order
    terms = [length(t.operators) == order for t in HT.terms]
    HT_order = SASQ.Expression(HT[terms])
    ex = P * HT_order |> act_on_ket |> simplify_heavy
    return ex
end

# Electronic cluster operator to n'th order
# T = 1/N! ∑_μ t_μ τ_μ
τ(n) = prod(E(2i-1,2i) for i = 1:n) * occupied(2:2:2n...) * virtual(1:2:2n...)
Tn(n) = 1 // factorial(n) * ∑(psym_tensor("t", 1:2n...) * τ(n), 1:2n)

# Define Hamiltonian in terms of F and g
Φ = 1//2 * ∑(psym_tensor("g", 1,2,3,4) * e(1,2,3,4), 1:4) +
    ∑(-2 * psym_tensor("g", 1,2,3,3) * E(1,2) * occupied(3), 1:3) +
    ∑( 1 * psym_tensor("g", 1,3,3,2) * E(1,2) * occupied(3), 1:3)

F = ∑(real_tensor("F", 1, 2) * E(1, 2), 1:2)

H = F + Φ

trans = translate(VirtualOrbital => 1:2:10, OccupiedOrbital => 2:2:10)

# Solve CCSDT equations (~2 min / 8 threads)
# Note: external symmetries ai<->bj<->ck is not found
T = Tn(2) + Tn(3)
Ecorr = cc_ket(H, T, 1, 0)
omega_ai = cc_ket(H, T, 2, 1)
omega_aibj = cc_ket(H, T, 2, 2)
omega_aibjck = cc_ket(H, T, 2, 3)

println("CC Energy:\n", (Ecorr, trans))
println("Ω_ai:\n", (omega_ai, trans))
println("Ω_aibj:\n", (omega_aibj, trans))
println("Ω_aibjck:\n", (omega_aibjck, trans))