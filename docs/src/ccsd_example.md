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

Here we see a very common pattern; the 2 * coulomb - 1 * exchange pattern where
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

### Omega Equations

#### Biorthogonal Bra

For the singles equations we want to compute quantity

``
\Omega_i^a
=
\bar{\left\langle_i^a\right|} \bar H |\text{HF}\rangle
``

where the bra ``\bar{\left\langle_i^a\right|}`` is the biorthogonal bra state
for the singly excited determinant

``\left|_i^a\right\rangle = E_{ai} |\text{HF}\rangle``

defined such that

``\left\langle\bar{_i^a}|_j^b\right\rangle = \delta_{ij} \delta_{ab}``

For a singly excited determinant we can explicitly express this biorthogonal bra
simply as

``\bar{\left\langle_i^a\right|} = \frac12 \langle\text{HF}| E_{ia}``

As a quick example we can show this

```@repl 1
bra = 1//2 * E(2, 1) * occupied(2) * virtual(1)
ket = E(3, 4) * occupied(4) * virtual(3)
hf_expectation_value(bra * ket)
```

However for higher excitations the biorthogonal bra states can not be as easily
expressed explicitly. For doubly excited determinants we want a biorthornormal
bra state ``\bar{\left\langle_{ij}^{ab}\right|}`` such that

``
\left\langle\bar{_{ij}^{ab}}|_{kl}^{cd}\right\rangle
=
P_{ai,bj} \delta_{ik} \delta_{jl} \delta_{ac} \delta_{bd}
=
\delta_{ik} \delta_{jl} \delta_{ac} \delta_{bd} +
\delta_{jk} \delta_{il} \delta_{bc} \delta_{ad}
``

with the permutation operator ``P_{ai,bj}`` generating the equivalent
permutations due to the particle symmetry in the doubly excited determinant

``
\left|_{ij}^{ab}\right\rangle =
\left|_{ji}^{ba}\right\rangle =
E_{ai} E_{bj} |\text{HF}\rangle =
E_{bj} E_{ai} |\text{HF}\rangle
``

We can write the biorthogonal bra explicitly as

``
\bar{\left\langle_{ij}^{ab}\right|}
=
\frac13 \left\langle_{ij}^{ab}\right| +
\frac16 \left\langle_{ji}^{ab}\right|
=
\langle\text{HF}| \left(
    \frac13 E_{jb} E_{ia} +
    \frac16 E_{ib} E_{ja}
\right)
``

Which we can confirm using the code

```@repl 1
bra = (1//3 * E(4, 3) * E(2, 1) + 1//6 * E(2, 3) * E(4, 1)) *
    occupied(2, 4) * virtual(1, 3)
ket = E(5, 6) * E(7, 8) * occupied(6, 8) * virtual(5, 7)
hf_expectation_value(bra * ket)
```

!!! note
    For higher order excited determinants such as triples, quadruples, etc. it is,
    however, not possible to explicitly express a biorthogonal state that has
    the desired properties. This is related to the fact that the expression is zero

    ``P^{abc} E_{ai} E_{bj} E_{ck} = 0``

    It is, however, still useful to be able to derive the expressions you would get
    if you project on a biorthogonal basis anyway, even if it is not expressible
    as a linear combination of excited determinants. For the triples for
    example, we can define two distinct expressions

    ``
    \Omega_{ijk}^{abc}
    =
    \left\langle\bar{_{ijk}^{abc}}\right| \bar H |\text{HF}\rangle
    ``

    ``
    \tilde \Omega_{ijk}^{abc}
    =
    \left\langle_{ijk}^{abc}\right| \bar H |\text{HF}\rangle
    ``

    the latter of which is the quantity that "actually exists" as a projection
    onto a set of bra states which are expressible as excited determinants and
    represent the equations which are to be solved to obtain the coupled cluster
    wave function. The two ``\Omega`` expressions are related by a
    non invertible linear transformation

    ``
    \tilde \Omega_{ijk}^{abc}
    =
    \sum_{a'i'b'j'c'k'}{
        S_{ijki'j'k'}^{abca'b'c'}
        \Omega_{i'j'k'}^{a'b'c'}
    }
    ``

    We can then see that if one finds amplitudes such that the ``\Omega``
    found by using the biorthogonal states they will also make the
    ``\tilde \Omega`` be zero, which is the condition we want to solve for.
    This is nice as we can safely "pretend" there exists a biorthogonal basis
    which produces much simpler expressions and are therefore more efficient to
    derive and numerically compute.

Since we do not want to (and often can not, see note above)
make an explicit expression for a biorthogonal bra we provide the function
[`project_biorthogonal`](@ref) which lets us insert Kronecker deltas according
to the definition of the biorthogonality in addition to the
[`symmetrize`](@ref) function to expand the permutations required by the
biorthogonality. We can reproduce the same behaviour as for the biorthogonal
doubles from above as

```@repl 1
ket = E(5, 6) * E(7, 8) * occupied(6, 8) * virtual(5, 7)
project_biorthogonal(ket, E(1, 2) * E(3, 4))
symmetrize(ans, make_permutation_mappings([(1, 2), (3, 4)]))
```

Here we supplied the "template" ket `E(1, 2) * E(3, 4)` which is that we want
to project on the biorthogonal bra of.

#### Singles Equations

We can use the [`project_biorthogonal`](@ref) function to project `Hbar_ket`
on the singles biorthogonal bra and get an expression for ``\Omega_i^a``

```@repl 1
omega_ai = project_biorthogonal(Hbar_ket, E(1, 2))
```

Here we also see the 2 * coulomb - 1 * exchange pattern show up again, so we can
simplify by.

```@repl 1
omega_ai = look_for_tensor_replacements(omega_ai,
    make_exchange_transformer("t", "u"))
```

Which produces a nice simplified expression for the singes part of the CCSD
equations.

#### Doubles Equations

Here we want to derive the expression for the doubles ``\Omega``

``
\Omega_{ij}^{ab}
=
\bar{\left\langle_{ij}^{ab}\right|} \bar H |\text{HF}\rangle
``

Just like for the singles equations we project on a biorthogonal bra, though
now using a doubly excited ket as the template, and we need to symmetrize
to include the pertmutations caused by the particle symmetry of the biorthogonal
bra.

```@repl 1
project_biorthogonal(Hbar_ket, E(1, 2) * E(3, 4));
omega_aibj = simplify_heavy(
    symmetrize(ans, make_permutation_mappings([(1, 2), (3, 4)])));
omega_aibj = look_for_tensor_replacements(omega_aibj,
    make_exchange_transformer("t", "u"))
```

Here we needed to use the (sometimes) expensive [`simplify_heavy`](@ref)
function to fully simplify, as well as recognizing the
2 * coloumb - 1 * exchange pattern for the T2 amplitudes.
Now we have a correct and rather nice expression for the doubles equations,
however, because of the particle symmetry of the bra that we introduced with
the call to [`symmetrize`](@ref) some terms are pairwise equal after permuting
the indices (a, i) <-> (b, j).
