# Quantum Electrodynamics CCSD Example

## Hamiltonian

In the QED-CC method for a single photon mode the Hamiltonian is of the form

``
H = H_e + \sum_{pq}{d_{pq} E_{pq} (b + b^\dagger)} +
\omega b^\dagger b
``

with ``H_e`` being the electronic Hamiltonian with adjusted integrals due to
the dipole self-energy, and ``b`` and ``b^\dagger`` are boson annihilation and
creation operators respectively for the cavity mode.

We define the electronic hamiltonian as before

```@repl 1
using SpinAdaptedSecondQuantization
h = ∑((
        real_tensor("F", 1, 2) +
        ∑((-2 * psym_tensor("g", 1, 2, 3, 3) +
            psym_tensor("g", 1, 3, 3, 2)) * occupied(3), [3])
    ) * E(1, 2) * electron(1, 2), 1:2);
g = 1//2 * simplify(
    ∑(psym_tensor("g", 1:4...) * e(1:4...) * electron(1:4...), 1:4)
);
He = h + g
```

Then the bilinear and pure photon terms

```@repl 1
b = boson()
H_bilinear = ∑(
    real_tensor("d", 1, 2) * E(1, 2) * (b' + b) * electron(1, 2), 1:2)
H_photon = real_tensor("ω") * b'b
H = He + H_bilinear + H_photon
```

## Cluster Operator

For the QED-CCSD method the cluster operator is given by

``
T = T1 + T2 + S0 + S1 + S2
``

with the separate parts defined as

```@repl 1
T2 = 1//2 * ∑(psym_tensor("t", 1:4...) * E(1, 2) * E(3, 4) *
occupied(2, 4) * virtual(1, 3), 1:4)
S0 = real_tensor("s0") * b'
S1 = ∑(real_tensor("s", 1, 2) * E(1, 2) * b' * occupied(2) * virtual(1), 1:2)
S2 = 1//2 * ∑(psym_tensor("s", 1:4...) * E(1, 2) * E(3, 4) * b' *
occupied(2, 4) * virtual(1, 3), 1:4)
```

As for electronic CC, we disregard `T1` as this is included in the integrals.

## Energy

```@repl 1
Hbar_ccsd = simplify(bch(He, T2, 4));
Hbar_ket_ccsd = simplify(act_on_ket(Hbar_ccsd, 2));
E_ccsd = simplify_heavy(act_on_bra(Hbar_ket_ccsd))
T = [T2, S0, S1, S2];
Hbar = simplify(bch(H, T, 4));
Hbar_ket = simplify(act_on_ket(Hbar, 3)); # Keeping 3 operators (E_ai E_bj b')
E0 = simplify_heavy(act_on_bra(Hbar_ket))
ΔE_QED_CCSD = E0 - E_ccsd # Additional terms from ccsd
```

## Omega

### Singles

```@repl 1
omega_ai_ccsd = project_biorthogonal(Hbar_ket_ccsd, E(1, 2));
omega_0_ai = project_biorthogonal(Hbar_ket, E(1, 2));
omega_0_ai - omega_ai_ccsd;
look_for_tensor_replacements(ans, make_exchange_transformer("t", "u"));
look_for_tensor_replacements(ans, make_exchange_transformer("s", "v"));
Δomega_ai = ans
```

### Doubles

```@repl 1
omega_aibj_ccsd = project_biorthogonal(Hbar_ket_ccsd, E(1, 2) * E(3, 4));
omega_0_aibj = project_biorthogonal(Hbar_ket, E(1, 2) * E(3, 4));
omega_0_aibj - omega_aibj_ccsd;
symmetrize(ans, make_permutation_mappings([(1, 2), (3, 4)]));
simplify_heavy(ans);
look_for_tensor_replacements(ans, make_exchange_transformer("t", "u"));
Δomega_aibj_r, Δomega_aibj_ss, Δomega_aibj_ns =
    desymmetrize(ans, make_permutation_mappings([(1, 2), (3, 4)]));
Δomega_aibj_r
Δomega_aibj_ss
Δomega_aibj_ns
```

### Photon

```@repl 1
project_biorthogonal(Hbar_ket, b');
look_for_tensor_replacements(ans, make_exchange_transformer("s", "v"));
omega_1 = ans
```

### Singles Photon

```@repl 1
project_biorthogonal(Hbar_ket, E(1, 2) * b');
simplify_heavy(ans);
look_for_tensor_replacements(ans, make_exchange_transformer("t", "u"));
look_for_tensor_replacements(ans, make_exchange_transformer("s", "v"));
look_for_tensor_replacements(ans, make_exchange_transformer("g", "L"));
omega_1_ai = ans
```

### Doubles Photon

```@repl 1
project_biorthogonal(Hbar_ket, E(1, 2) * E(3, 4) * b');
symmetrize(ans, make_permutation_mappings([(1, 2), (3, 4)]));
simplify_heavy(ans);
look_for_tensor_replacements(ans, make_exchange_transformer("t", "u"));
look_for_tensor_replacements(ans, make_exchange_transformer("s", "v"));
look_for_tensor_replacements(ans, make_exchange_transformer("g", "L"));
omega_1_aibj_r, omega_1_aibj_ss, omega_1_aibj_ns =
    desymmetrize(ans, make_permutation_mappings([(1, 2), (3, 4)]));
omega_1_aibj_r
omega_1_aibj_ss
omega_1_aibj_ns
```
