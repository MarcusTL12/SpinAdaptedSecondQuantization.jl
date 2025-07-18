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

!!! note
    Usually one wants to evaluate projections using the *biorthonormal* bra
    which for doubles is defined as

    ``
    \tilde{\left\langle_{ij}^{ab}\right|} =
    \frac{1}{1 + \delta_{ai,bj}}
    \bar{\left\langle_{ij}^{ab}\right|}
    ``

    which removes the factor of 2 on the diagonal elements where `ai = bj`.
    This is usually done by numerically scaling the diagonal of the output by
    0.5 as the expressions look nicer using the biorthogonal bra. Note that for
    higher order excitations there are more "diagonals" to think about.

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
omega_aibj_u = look_for_tensor_replacements(omega_aibj,
    make_exchange_transformer("t", "u"))
```

Here we needed to use the (sometimes) expensive [`simplify_heavy`](@ref)
function to fully simplify, as well as recognizing the
2 * coulomb - 1 * exchange pattern for the T2 amplitudes.
Now we have a correct and rather nice expression for the doubles equations,
however, because of the particle symmetry of the bra that we introduced with
the call to [`symmetrize`](@ref) some terms are equal after permuting
the indices (a, i) <-> (b, j). Computing all of these redundant terms leads to
a lot of unnecessary terms to code and compute, so it is much more efficient
to only code the non-redundant ones, then symmetrize the final result
numerically. To help with this, we provide the [`desymmetrize`](@ref) function
which sort of acts like the reverse of [`symmetrize`](@ref). We can call this
function on our ``\Omega_{ij}^{ab}`` expression

```@repl 1
omega_aibj_r, omega_aibj_ss, omega_aibj_ns = desymmetrize(
    omega_aibj_u, make_permutation_mappings([(1, 2), (3, 4)])
);
omega_aibj_r
omega_aibj_ss
omega_aibj_ns
```

The [`desymmetrize`](@ref) function returns three expressions. The first
`omega_aibj_r` contains the terms it found redundant permutations of elsewhere
in the original expression. These needs to be symmetrized to obtain the correct
expression. The second term `omega_aibj_ss` contains the "self-symmetric" terms.
These are the terms which themselves carry the symmetry requested and can be
evaluated as-is. The third and last term `omega_aibj_ns` are the "non-symmetric"
terms. Interestingly, we have here gotten two terms in the non-symmetric
expression, even though the result should be symmetric. Under closer inspection
one can see that the two "non-symmetric" terms are indeed symmetric if one
re-expands the ``u_{aidk}`` tensor in terms of ``t``. To avoid having
non-symmetric terms we can instead directly desymmetrize the omega before
looking for coulomb - exchange symmetry.

```@repl 1
omega_aibj_r, omega_aibj_ss, omega_aibj_ns = desymmetrize(
    omega_aibj, make_permutation_mappings([(1, 2), (3, 4)])
);
look_for_tensor_replacements(omega_aibj_r, make_exchange_transformer("g", "L"))
look_for_tensor_replacements(omega_aibj_ss, make_exchange_transformer("g", "L"))
omega_aibj_ns
```

where we see that the non-symmetric part is indeed zero, but we have lost the
ability to remove all prefactors. Counting total terms we can see that both
strategies give 15 terms in total, so the cost should not be much different, but
it can be nice to keep in mind. In any case to evaluate the full doubles omega
one needs to symmetrize the redundant expression, and add the direct evaluation
of the self-symmetric and apparent non-symmetric expression which we can write
as

``
\Omega_{ij}^{ab} =
(\Omega_{\text{ss}})_{ij}^{ab} +
(\Omega_{\text{ns}})_{ij}^{ab} +
P_{ij}^{ab} (\Omega_{\text{s}})_{ij}^{ab}
``

## Jacobian Transformation

The coupled cluster Jacobian is a central quantity for computing various
properties of the coupled cluster wave function, such as density matrices and
excitation energies. The Jacobian matrix ``A_{\mu\nu}`` is given by

``
A_{\mu\nu} = \langle\mu| e^{-T} [H, \tau_{\nu}] e^T |\text{HF}\rangle
``

though one is rarely interested in the full matrix, but rather to compute
eigenvalues and solve linear systems. The matrix is usually huge and very
structures (though not very sparse), so it is way more efficient to code a
direct linear transformation. The Jacobian is not symmetric so one needs to
derive both the right and left transformations given by

``
\rho_\mu = \sum_{\nu}{A_{\mu\nu} c_\nu}
``

``
\sigma_\nu = \sum_{\mu}{b_\mu A_{\mu\nu}}
``

We can quite nicely derive expressions for both of these transformations using
the code.

### Right transformation

The transformed vector ``\rho_\mu`` can be split into singles and doubles
block like

``
\rho_{ai} = \sum_{\nu}{A_{ai,\nu} c_\nu}
``

``
\rho_{aibj} = \sum_{\nu}{A_{aibj,\nu} c_\nu}
``

#### Singles block

We then split the sum into blocks

``
\rho_{ai} = \sum_{\nu}{A_{ai,\nu} c_\nu}
=
\sum_{ck}{A_{ai,ck} c_{ck}} +
\sum_{ck \geq dl}{A_{ai,ckdl} c_{ckdl}}
``

The singles sum is okay, but the triangular doubles sum ``ck \geq dl`` is not
something we can express using the sums in the code, so we have to square it up.

``
\sum_{ck \geq dl}{A_{ai,ckdl} c_{ckdl}}
=
\frac12 \sum_{ckdl}{(1 + \delta_{ck,dl}) A_{ai,ckdl} c_{ckdl}}
=
\frac12 \sum_{ckdl}{A_{ai,ckdl} \tilde c_{ckdl}}
``

Here we have to add a factor of 2 on the diagonal of ``c_{ckdl}`` by defining
the adjusted tensor

``\tilde c_{ckdl} = (1 + \delta_{ck,dl}) c_{ckdl}``

though we will not be using the ``\tilde c`` name in the derivation, but rather
remembering that we have to numerically scale the diagonal of ``c_{ckdl}`` by
2 when coding the final expression.

We start by computing the projection

``
\sum_{ck}{e^{-T} [H, E_{ck}] e^T |\text{HF}\rangle c_{ck}}
``

```@repl 1
simplify(commutator(H, E(1, 2) * occupied(2) * virtual(1)));
simplify(bch(ans, T2, 4));
proj_ai = simplify(act_on_ket(ans, 2));
proj_ck = simplify(∑(ans * real_tensor("c", 1, 2), [1, 2]));
```

This projection will also be useful when doing the doubles block.
We also keep the `proj_ai` for the left transformation later.
To get the singles block out we project on a biorthogonal singles bra

```@repl 1
ρ_ai_singles = project_biorthogonal(proj_ck, E(1, 2));
ρ_ai_singles = look_for_tensor_replacements(ρ_ai_singles,
    make_exchange_transformer("t", "u"));
ρ_ai_singles = look_for_tensor_replacements(ρ_ai_singles,
    make_exchange_transformer("g", "L"));
ρ_ai_singles
```

where we get a very consise expression for the singles-singles transformation.
Next we compute the projection

``
\frac12 \sum_{ckdl}{e^{-T} [H, E_{ck} E_{dl}] e^T |\text{HF}\rangle c_{ckdl}}
``

```@repl 1
simplify(commutator(H, E(1, 2) * E(3, 4) * occupied(2, 4) * virtual(1, 3)));
simplify(bch(ans, T2, 4));
proj_aibj = simplify(act_on_ket(ans, 2));
proj_ckdl = simplify(1//2 * ∑(ans * psym_tensor("c", 1:4...), 1:4));
```

This will also be useful for the doubles block later and we also keep the
proj_aibj for the left transformation.
Now we project on the singles bra to get

```@repl 1
ρ_ai_doubles = project_biorthogonal(proj_ckdl, E(1, 2));
ρ_ai_doubles = look_for_tensor_replacements(ρ_ai_doubles,
    make_exchange_transformer("g", "L"));
ρ_ai_doubles = look_for_tensor_replacements(ρ_ai_doubles,
    make_exchange_transformer("c", "cu"));
ρ_ai_doubles
```

Here we introduced another tensor pattern

``cu_{aibj} = 2 c_{aibj} - c_{ajbi}``

We can now combine the two parts of the ``\rho_{ai}`` block

```@repl 1
ρ_ai = ρ_ai_singles + ρ_ai_doubles
```

#### Doubles block

We again split the sum into singles and doubles, squareup and include the
factor of 2 on the diagonal

``
\rho_{aibj} = \sum_{\nu}{A_{aibj,\nu} c_\nu}
=
\sum_{ck}{A_{aibj,ck} c_{ck}} +
\frac12 \sum_{ckdl}{A_{aibj,ckdl} \tilde c_{ckdl}}
``

Here we can reuse the right projection from above, projecting on the
biorthogonal doubles bra instead.

First the singles sum

```@repl 1
project_biorthogonal(proj_ck, E(1, 2) * E(3, 4));
symmetrize(ans, make_permutation_mappings([(1, 2), (3, 4)]));
simplify_heavy(ans);
look_for_tensor_replacements(ans, make_exchange_transformer("t", "u"));
look_for_tensor_replacements(ans, make_exchange_transformer("g", "L"));
ρ_aibj_singles_r, ρ_aibj_singles_ss, ρ_aibj_singles_ns =
    desymmetrize(ans, make_permutation_mappings([(1, 2), (3, 4)]));
ρ_aibj_singles_r
ρ_aibj_singles_ss
ρ_aibj_singles_ns
```

Then the doubles sum


```@repl 1
project_biorthogonal(proj_ckdl, E(1, 2) * E(3, 4));
symmetrize(ans, make_permutation_mappings([(1, 2), (3, 4)]));
simplify_heavy(ans);
look_for_tensor_replacements(ans, make_exchange_transformer("t", "u"));
look_for_tensor_replacements(ans, make_exchange_transformer("g", "L"));
ρ_aibj_doubles_r, ρ_aibj_doubles_ss, ρ_aibj_doubles_ns =
    desymmetrize(ans, make_permutation_mappings([(1, 2), (3, 4)]));
ρ_aibj_doubles_r
ρ_aibj_doubles_ss
ρ_aibj_doubles_ns
```

Combining singles and doubles

```@repl 1
ρ_aibj_r = ρ_aibj_singles_r + ρ_aibj_doubles_r
ρ_aibj_ss = ρ_aibj_doubles_ss
```

which concludes the jacobian right transformation.

### Left transformation

The procedure for the left transformation is very similar to the right
transformation with some key differences. We can still split the transformed
vector into singles and doubles blocks

The transformed vector ``\rho_\mu`` can be split into singles and doubles
block like

``
\sigma_{ai} = \sum_{\mu}{b_\mu A_{\mu,ai}}
``

``
\sigma_{aibj} = \sum_{\mu}{b_\mu A_{\mu,aibj}}
``

#### Singles block

We again further split the sum into singles and doubles parts

``
\sigma_{ai} = \sum_{\mu}{b_\mu A_{\mu,ai}}
=
\sum_{ck}{b_{ck} A_{ck,ai}} +
\frac12 \sum_{ckdl}{(1 + δ_{ck,dl}) b_{ckdl} A_{ckdl,ai}}
``

though now we get a nice cancellation due to the fact that the jacobian
is defined in terms of the *biorthonormal* bra, so we can define an adjusted
jacobian

``\tilde A_{ckdl,\mu} = \frac{1}{1 + δ_{ck,dl}} A_{ckdl,\mu}``

defined in terms of the biorthogonal bra. The factors of 2 on the diagonal
thus cancel and we are left with the much nicer expression

``
\sigma_{ai}
=
\sum_{ck}{b_{ck} A_{ck,ai}} +
\frac12 \sum_{ckdl}{b_{ckdl} \tilde A_{ckdl,ai}}
``

where we do not have to scale either the input vector ``b`` or the output
vector ``\sigma`` on the diagonals, contrary to the right transformation.

We can reuse the projection from above and can directly compute the singles
projection

```@repl 1
simplify(
    ∑(project_biorthogonal(proj_ai, E(5, 6)) * real_tensor("b", 5, 6), 5:6)
);
look_for_tensor_replacements(ans, make_exchange_transformer("t", "u"));
look_for_tensor_replacements(ans, make_exchange_transformer("g", "L"));
σ_ai_singles = ans
```

and similrarly for the doubles

```@repl 1
project_biorthogonal(proj_ai, E(5, 6) * E(7, 8)) * psym_tensor("b", 5:8...);
simplify_heavy(∑(ans, 5:8)); # No need to symmetrize because of sum
look_for_tensor_replacements(ans, make_exchange_transformer("t", "u"));
look_for_tensor_replacements(ans, make_exchange_transformer("g", "L"));
σ_ai_doubles = ans
```

!!! note
    Here we omited the call to `symmetrize` after `project_biorthogonal` because
    simplify can reorder the indices in the sum and achieve the same. This
    misses a factor of 2 which is compensated for by omiting the `1//2` in front
    of the sum.

And combined

```@repl 1
σ_ai = σ_ai_singles + σ_ai_doubles
```

#### Doubles block

The same cancellations hold as above, so the procedure is almost identical, but
using the other ket projection.

Singles:

```@repl 1
simplify(
    ∑(project_biorthogonal(proj_aibj, E(5, 6)) * real_tensor("b", 5, 6), 5:6)
);
look_for_tensor_replacements(ans, make_exchange_transformer("g", "L"));
σ_aibj_singles_r, σ_aibj_singles_ss, σ_aibj_singles_ns =
    desymmetrize(ans, make_permutation_mappings([(1, 2), (3, 4)]));
σ_aibj_singles_r
σ_aibj_singles_ss
σ_aibj_singles_ns
```

Doubles:

```@repl 1
project_biorthogonal(proj_aibj, E(5, 6) * E(7, 8)) * psym_tensor("b", 5:8...);
simplify_heavy(∑(ans, 5:8)); # No need to symmetrize because of sum
look_for_tensor_replacements(ans, make_exchange_transformer("t", "u"));
look_for_tensor_replacements(ans, make_exchange_transformer("g", "L"));
σ_aibj_doubles_r, σ_aibj_doubles_ss, σ_aibj_doubles_ns =
    desymmetrize(ans, make_permutation_mappings([(1, 2), (3, 4)]));
σ_aibj_doubles_r
σ_aibj_doubles_ss
σ_aibj_doubles_ns
```

Combined:

```@repl 1
σ_aibj_r = σ_aibj_singles_r + σ_aibj_doubles_r
σ_aibj_ss = σ_aibj_doubles_ss
```

Which concludes the left transformation.

## One-electron density matrices

One frequently wants the one-electron density matrix to get properties like
dipole, quadrupole, magnetic moments, etc. both for a state and
transition moments. In coupled cluster theory we have asymmetric density
matrices which we can express as

``D_{pq} = \langle L| E_{pq} |R\rangle``

where the left and right states are given in general by

``
\langle L|
=
\left(
    L_0 \langle\text{HF}| +
    \sum_{\mu}{L_\mu \langle\mu|}
\right) e^{-T}
``

``
|R\rangle
=
e^T \left(
    |\text{HF}\rangle R_0
    +
    \sum_{\nu}{|\nu\rangle R_\nu}
\right)
``

The left and right coefficient vectors can represent either ground excited
states giving us the possibility to derive general expressions for ground
and excited state as well as transition densities.

Firstly we will need the similarity transformed operator for all parts

```@repl 1
bch_Epq = bch(E(1, 2) * electron(1, 2), T2, 4)
```

### Ref-Ref contributions

The terms arising from the ``L_0`` and ``R_0`` (reference-reference)
are quite straight forward to derive:

``
D_{pq}^\text{ref-ref}
=
L_0 \langle\text{HF}| e^{-T} E_{pq} e^T |\text{HF}\rangle R_0
``

```@repl 1
real_tensor("L0") * real_tensor("R0") * bch_Epq;
hf_expectation_value(ans)
```

where we see that only the HF density survives.

!!! note
    Since we are not including T1 in the cluster operator we will initially
    miss many terms in the density coming from them. This is fine if the
    density is contracted with integrals that are T1 transformed, or one can
    include the T1 terms into the density by a similarity transform numerically
    which makes the expressions a bit nicer.

### ``\mu``-Ref contributions

The terms arising from ``L_\mu`` and ``R_0`` are the only additional terms
needed for the ground state density as ``R_\nu = 0`` for the ground state.
We will omit the ``R_0`` in the derivation to reduce clutter.
We can split the term into singles and doubles parts

``
D_{pq}^{\mu\text{-ref}}
=
\sum_{\mu}{L_\mu \langle\mu|} e^{-T} E_{pq} e^T |\text{HF}\rangle
\\=
\sum_{ai}{L_{ai} \tilde{\langle_i^a|}} e^{-T} E_{pq} e^T |\text{HF}\rangle +
\sum_{ai \geq bj}{L_{aibj} \tilde{\langle_{ij}^{ab}|}}
e^{-T} E_{pq} e^T |\text{HF}\rangle
``

in both cases we want the projection

``e^{-T} E_{pq} e^T |\text{HF}\rangle``

```@repl 1
ref_ket = simplify_heavy(act_on_ket(bch_Epq, 2));
```

Singles:

```@repl 1
∑(project_biorthogonal(ref_ket, E(3, 4)) * real_tensor("L", 3, 4), 3:4);
simplify(ans);
look_for_tensor_replacements(ans, make_exchange_transformer("t", "u"))
```

For doubles we have the same cancellation as for the Jacobian left
transformation

``
\sum_{ai \geq bj}{L_{aibj} \tilde{\langle_{ij}^{ab}|}}
e^{-T} E_{pq} e^T |\text{HF}\rangle
=
\frac12 \sum_{aibj}{L_{aibj} \bar{\langle_{ij}^{ab}|}}
e^{-T} E_{pq} e^T |\text{HF}\rangle
``

```@repl 1
project_biorthogonal(ref_ket, E(3, 4) * E(5, 6));
∑(ans * psym_tensor("L", 3:6...), 3:6);
D_mu_ref = simplify_heavy(ans)
```

!!! note
    Since the two terms in this expression have external indices with different
    constraints, it can be a bit confusing to look at. Since these will be
    coded as separate blocks of the density matrix, it is usually advisable
    to extract the various blocks into their own expressions

    ```@repl 1
    D_ij = D_mu_ref * occupied(1, 2)
    D_ab = D_mu_ref * virtual(1, 2)
    ```

    Alternatively one can disable translation of external indices to see
    the actual indices and constraints
    ```@repl 1
    disable_external_index_translation()
    D_mu_ref # Here we see index ₁ and ₂ as we specified in E(1, 2)
    enable_external_index_translation()
    ```

### Ref-``\nu`` contributions

``
D_{pq}^{\text{ref-}\nu}
=
L_0 \sum_{\nu}{\langle\text{HF}| e^{-T} E_{pq} e^T |\nu\rangle R_\nu}
\\=
L_0 \sum_{ai}{\langle\text{HF}| e^{-T} E_{pq} e^T |_i^a\rangle R_{ai}} +
L_0 \sum_{ai \geq bj}{
    \langle\text{HF}| e^{-T} E_{pq} e^T |_{ij}^{ab}\rangle R_{aibj}
}
``

Singles:

```@repl 1
∑(bch_Epq * E(3, 4) * occupied(4) * virtual(3) * real_tensor("R", 3, 4), 3:4);
singles_ket = simplify(act_on_ket(ans, 2));
simplify_heavy(act_on_bra(singles_ket))
```

For doubles we scale the ``R_{aibj}`` by 2 on the diagonal

``
\sum_{ai \geq bj}{
    \langle\text{HF}| e^{-T} E_{pq} e^T |_{ij}^{ab}\rangle R_{aibj}
}
=
\frac12 \sum_{aibj}{
    (1 + \delta_{ai,bj})
    \langle\text{HF}| e^{-T} E_{pq} e^T |_{ij}^{ab}\rangle R_{aibj}
}
=
\frac12 \sum_{aibj}{
    \langle\text{HF}| e^{-T} E_{pq} e^T |_{ij}^{ab}\rangle \tilde R_{aibj}
}
``

```@repl 1
bch_Epq * E(3, 4) * E(5, 6) * occupied(4, 6) * virtual(3, 5);
simplify(1//2 * ∑(ans * psym_tensor("R", 3:6...), 3:6));
doubles_ket = simplify(act_on_ket(ans, 2));
simplify_heavy(act_on_bra(doubles_ket))
```

### ``\mu``-``\nu`` contributions

``
D_{pq}^{\mu\text{-}\nu}
=
\sum_{\mu\nu}{
    L_\mu \langle\mu| e^{-T} E_{pq} e^T |\nu\rangle R_\nu
}
\\=
\sum_{aick}{
    L_{ai} \tilde{\langle_i^a|} e^{-T} E_{pq} e^T |_k^c\rangle R_{ck}
} \\+
\sum_{aibjck}{
    L_{aibj} \tilde{\langle_{ij}^{ab}|} e^{-T} E_{pq} e^T |_k^c\rangle R_{ck}
} \\+
\sum_{aickdl}{
    L_{ai} \tilde{\langle_i^a|} e^{-T} E_{pq} e^T |_{kl}^{cd}\rangle R_{ckdl}
} \\+
\sum_{aibjckdl}{
    L_{aibj} \tilde{\langle_{ij}^{ab}|} e^{-T} E_{pq} e^T
    |_{kl}^{cd}\rangle R_{ckdl}
}
``

Singles-Singles:

```@repl 1
project_biorthogonal(singles_ket, E(3, 4));
simplify(∑(ans * real_tensor("L", 3, 4), 3:4));
D_ss = ans
```

Doubles-Singles:

```@repl 1
project_biorthogonal(singles_ket, E(3, 4) * E(5, 6));
simplify_heavy(∑(ans * psym_tensor("L", 3:6...), 3:6));
look_for_tensor_replacements(ans, make_exchange_transformer("t", "u"));
D_ds = ans
```

Singles-Doubles:

```@repl 1
project_biorthogonal(doubles_ket, E(3, 4));
simplify(∑(ans * real_tensor("L", 3, 4), 3:4));
look_for_tensor_replacements(ans, make_exchange_transformer("R", "Ru"));
D_sd = ans
```

Doubles-Doubles:

```@repl 1
project_biorthogonal(doubles_ket, E(3, 4) * E(5, 6));
simplify_heavy(∑(ans * psym_tensor("L", 3:6...), 3:6));
D_dd = ans
```

Combining and separating into blocks:

```@repl 1
D = D_ss + D_sd + D_ds + D_dd
D_ij = D * occupied(1, 2)
D_ia = D * occupied(1) * virtual(2)
D_ai = D * occupied(2) * virtual(1)
D_ab = D * virtual(1, 2)
```
