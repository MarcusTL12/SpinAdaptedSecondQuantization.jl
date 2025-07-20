
# Deriving Box 13.2 from "Molecular Electronic Structure Theory"

Here we show how to reproduce box 13.2 from the book
[Molecular Electronic Structure Theory](https://books.google.no/books?hl=en&lr=&id=APjLWFFxWkQC&oi=fnd&pg=PR21&dq=molecular+electronic+structure+theory&ots=Hxpe4uwYLg&sig=lgjszA3r6d9A1xVAxTXlVX5T_uY&redir_esc=y#v=onepage&q=molecular%20electronic%20structure%20theory&f=false)

Setup:

```@repl 1
using SpinAdaptedSecondQuantization
h = ∑((
        real_tensor("F", 1, 2) +
        ∑((-2 * psym_tensor("g", 1, 2, 3, 3) +
            psym_tensor("g", 1, 3, 3, 2)) * occupied(3), [3])
    ) * E(1, 2) * electron(1, 2), 1:2);
g = 1//2 * ∑(psym_tensor("g", 1:4...) * e(1:4...) * electron(1:4...), 1:4);
H = h + g
Eai(a, i) = E(a, i) * virtual(a) * occupied(i)
```

Entry 1:

```@repl 1
simplify_heavy(act_on_ket(H))
```

Since we do not have both the one electron matrix ``h_{pq}`` and the Fock matrix
``F_{pq}`` we do not perfectly reproduce the form of the HF energy from
box 13.2, but the rest of the box is reproduced nicely.

Entry 2:

```@repl 1
simplify_heavy(act_on_ket(commutator(H, Eai(1, 2))));
look_for_tensor_replacements(ans, make_exchange_transformer("g", "L"))
```

Entry 3:

```@repl 1
simplify_heavy(act_on_ket(commutator(H, Eai(1, 2), Eai(3, 4))));
look_for_tensor_replacements(ans, make_exchange_transformer("g", "L"));
r, ss, ns = desymmetrize(ans, make_permutation_mappings([(1, 2), (3, 4)]));
r # Terms that have P_ij^ab in front
ss # Symmetric terms
ns # Non symmetric terms
```

Entry 4:

```@repl 1
simplify_heavy(act_on_ket(commutator(H, Eai(1, 2), Eai(3, 4), Eai(5, 6))));
look_for_tensor_replacements(ans, make_exchange_transformer("g", "L"));
r, ss, ns = desymmetrize(ans,
    make_permutation_mappings([(1, 2), (3, 4), (5, 6)]));
r # Terms that have P_ijk^abc in front
ss # Symmetric terms
ns # Non symmetric terms
```

Entry 5:

```@repl 1
x = simplify_heavy(
    act_on_ket(commutator(H, Eai(1, 2), Eai(3, 4), Eai(5, 6), Eai(7, 8))));
r, ss, ns = desymmetrize(x,
    make_permutation_mappings([(1, 2), (3, 4), (5, 6)]));
r # The first part with P_ijk^abc
r, ss, ns = desymmetrize(x,
    make_permutation_mappings([(1, 2), (3, 4), (5, 6), (7, 8)]));
r # The second part with P_ijkl^abcd
```
