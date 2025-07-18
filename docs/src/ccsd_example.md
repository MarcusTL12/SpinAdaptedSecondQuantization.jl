# Coupled Cluster Singles and Doubles Example

Here is a detailed example of how to derive common expressions and equations for
the closed shell CCSD method. This will also serve as a way of showcasing the
various functionality of the package.

## Hamiltonian

For the hamiltonian we will use the familiar expression for the non-relativistic
electronic hamiltonian, however we will express it in terms of the fock matrix
instead of the normal one-electron matrix

``
F_{pq} = h_{pq} + \sum_{i}{\left(2 g_{pqii} - g_{piiq}\right)}
``

``
h_{pq} = F_{pq} - \sum_{i}{\left(2 g_{pqii} - g_{piiq}\right)}
``

This is because every expression except the HF energy looks simpler expressed
in terms of the fock matrix. In code we then define the hamiltonian as

```@repl 1
using SpinAdaptedSecondQuantization
h = ∑((
        real_tensor("F", 1, 2) +
        ∑((-2 * psym_tensor("g", 1, 2, 3, 3) +
            psym_tensor("g", 1, 3, 3, 2)) * occupied(3), [3])
    ) * E(1, 2) * electron(1, 2), 1:2)
g = 1//2 * simplify(
    ∑(psym_tensor("g", 1:4...) * e(1:4...) * electron(1:4...), 1:4)
)
H = h + g
```

We omit the nuclear repulsion term and we do not assume any symmetry in our
integrals apart from the particle exchange symmetry in the two electron
integrals. We do this since we can then use the T1 transformed integrals
allowing us to ignore the T1 operator when deriving equations.

We can derive the HF energy expression for the Hamiltonian in terms of the
Fock matrix as

```@repl 1
E_HF = simplify_heavy(hf_expectation_value(H))
```

## Cluster Operator

The cluster operator we define then as only the T2 operator as the T1 amplitudes
are included in the T1 transformed integrals.

```@repl 1
T2 = 1//2 * ∑(psym_tensor("t", 1:4...) * E(1, 2) * E(3, 4) *
occupied(2, 4) * virtual(1, 3), 1:4)
```

## Similrarity transformed Hamiltonian

A central operator in coupled cluster theory is the similarity transformed
Hamiltonian operator

``\bar H = e^{-T} H e^T``

This is usually computed using the Baker-Campbell-Hausdorff (BCH) expansion

``
e^{-B} A e^B
=
A + [A, B] +
\frac12 [[A, B], B] +
\frac{1}{3!} [[[A], B], B] +
...
``

A nice feature of this expansion for the electronic hamiltionian and the cluster
operator is that it truncates after the 4th order term (5 terms).
This is because the 5th order expression is zero

``
[[[[[H, E_{ai}], E_{bj}], E_{ck}], E_{dl}], E_{em}] = 0
``

This means that we can explicitly compute the full similarity transformed
Hamiltonian with no fear of missing any terms. We supply the function
`bch(A, B, n)` which computes the BCH expansion to nth order.

!!! note
    We also supply the functon `bch(A, [B1, B2, ...], n)` which
    is equivalent to `bch(A, B1 + B2 + ..., n)` when all of `B1, B2, ...`
    comute among each other, but is much faster as it avoids computing
    equivalent terms such as ``[[A, B_1], B_2]`` and ``[[A, B_2], B_1]``.
    This can be especially useful for higher order coupled cluster methods such
    as CCSDT and beyond, as well as exotic methods such as QED-CC and NEO-CC.

In code this looks like

```@repl 1
Hbar = simplify(bch(H, T2, 4)); # Avoid printing as this is quite long
```

## Ground State Energy Expression and Equations

The coupled cluster ground state energy and equations are obtained by projecting
the similarity constrained Schrödinger equation

``
\bar H |\text{HF}\rangle = E |\text{HF}\rangle
``

onto a set of bra states. The energy is obtained by projecting on the HF bra

``
E_0 = \langle\text{HF}| \bar H |\text{HF}\rangle
``

while the coupled cluster wave function is determined by projecting on a
set of excited bra states determined from the cluster operator

``
\langle\mu| \bar H |\text{HF}\rangle = 0
``

In both cases we require the quantity ``\bar H |\text{HF}\rangle`` which we
can explicitly compute by the function `act_on_ket(ex, [max_ops])`.

```@repl 1
Hbar_ket = simplify(act_on_ket(Hbar, 2)); # Avoid printing long expression
```

Here we choose to throw away any terms with more than two operators as in CCSD
we are never projecting on any bra that is more than doubly excited. This can
speed up the projection and further operations by a bit.

### Energy

To obtain the energy expression we project the Hbar_ket on HF from the left.
This can be nicely done using the `act_on_bra(ex, [max_ops])` function

```@repl 1
E0 = simplify_heavy(act_on_bra(Hbar_ket))
```

Since the interesting part here is that which differs from HF, we can get the
correlation energy as

```@repl 1
E_corr = E0 - E_HF
```

Here we see a very common pattern; the 2 * coulom - 1 * exchange pattern where
we have terms that differ by a prefactor of ``-\frac12`` and an exchange of
index 2 and 4 of a 4 index tensor, such as the two electron integrals or the
T2 amplitudes. We usually define separate tensor names for these patterns 
such as

``
L_{pqrs} = 2 g_{pqrs} - g_{psrq}
``

and

``
u_{aibj} = 2 t_{aibj} - t_{ajbi}
``

We provide a pair of functions to look for this pattern in expressions.
firstly and mainly [`look_for_tensor_replacements`](@ref) function
which is the actual engine that takes an expression and a tensor transformation
function and looks for terms that are equal up to this transformation.
We provide the utility function [`make_exchange_transformer`](@ref)
which produces a transformer that looks for the 2 * coulomb - 1 * exchange
pattern. We can show this in action on the correlation energy looking for both
the patterns mentioned above

```@repl 1
look_for_tensor_replacements(E_corr, make_exchange_transformer("g", "L"))
look_for_tensor_replacements(E_corr, make_exchange_transformer("t", "u"))
```
