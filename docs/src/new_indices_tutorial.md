
# Defining new index spaces

Defining your own index spaces can be a powerful tool when working with methods
where you are dealing with other types of indices than the default provided
spaces `GeneralOrbital`, `OccupiedOrbital` and `VirtualOrbital`.
A few examples could be
- Auxiliary index for integral decompositions (RI, Cholesky)
- Multiple photon modes for QED methods
- Subdividing `GeneralOrbital` to have different classes of fermions,
    for example:
    - Dividing into active and inactive orbitals for active space methods
    - Dividing into orbitals for different kinds of physical fermions such as
        protons and positrons (NEO-CCSD)

## RI/Cholesky index

The first example is adding a simple index space for use as an auxiliary basis
for integral decompositions. We define a new index space using the
[`new_space`](@ref) function

```@repl 1
using SpinAdaptedSecondQuantization
AuxiliaryIndex = new_space(:AuxiliaryIndex, "aux", "JKLMN")
```

The variable `AuxiliaryIndex` now contains the id for the new index space which
are ordered indices starting at `1`. This can also be obtained from the global
dictionary `index_ids`

```@repl 1
SASQ.index_ids[:AuxiliaryIndex]
```

The second parameter `"aux"` is the "short name" for the space, used when
printing terms without index translation. The last parameter `"JKLMN"` are the
names to give indices in the given space when translated.

!!! note
    Running `new_space` for the same key multiple times in the same session will
    not do anything, meaning it is completely safe to keep index space
    definitions in the top of a script that gets run multiple times.

We can now use our index space

```@repl 1
x = real_tensor("L", 1, 2, 3) * constrain(1 => AuxiliaryIndex) * electron(2, 3)
disable_external_index_translation()
x # Here the short name "aux" gets used
enable_external_index_translation()
```

It can be nice to define a convenience function to constrain many indices
to the new space

```@repl 1
auxiliary(indices...) = constrain(J => AuxiliaryIndex for J in indices)
real_tensor("S", 1, 2) * auxiliary(1, 2)
```

We can now, for example, define the two electron operator in terms of
auxiliary indices

```@repl 1
g = 1//2 * ∑(real_tensor("L", 5, 1, 2) * real_tensor("L", 5, 3, 4) *
    e(1, 2, 3, 4) * electron(1, 2, 3, 4) * auxiliary(5), 1:5)
simplify_heavy(hf_expectation_value(g))
```

## Multilevel indices

Say we want to divide our electronic indices into an active and an inactive
space. Then we can define the following spaces.

```@repl 1
ActiveOrbital = new_space(:ActiveOrbital, "a", "pqrstuv")
InactiveOrbital = new_space(:InactiveOrbital, "i", "pqrstuv")
```

Note that the same names have been given to these indices as the general
electronic indices, so they will be printed the same when translated. To be
able to distinguish them we can add a color to be printed for a certain index
space (see [`set_color`](@ref) for available colors)

```@repl 1
set_color(ActiveOrbital, :cyan)
set_color(InactiveOrbital, :red)
```

now if we use the indices we see that they get possible to tell apart

```@repl 1
∑(real_tensor("g", 1, 2, 3, 4) *
    electron(1, 2) * constrain(3 => ActiveOrbital, 4 => InactiveOrbital), 1:4)
```

However, we want these two new spaces to act in such a way that
```
ActiveOrbital ⊆ GeneralOrbital
InctiveOrbital ⊆ GeneralOrbital
ActiveOrbital ∪ InctiveOrbital == GeneralOrbital
ActiveOrbital ∩ InctiveOrbital == {}
```

To setup these relations we call the function [`add_space_sum`](@ref)

```@repl 1
add_space_sum(ActiveOrbital, InactiveOrbital, GeneralOrbital)
ActiveOrbital ⊆ GeneralOrbital
InactiveOrbital ⊆ GeneralOrbital
```

This for example allows the following simplification to occur

```@repl 1
∑(real_tensor("h", 1, 2) * E(1, 2) *
    electron(1) * constrain(2 => ActiveOrbital), 1:2)
ans + ∑(real_tensor("h", 1, 2) * E(1, 2) *
    electron(1) * constrain(2 => InactiveOrbital), 1:2)
simplify(ans)
```

Now we also want to subdivide the `OccupiedOrbital` and `OccupiedOrbital`
index spaces into active and inactive subspaces, so we define these spaces

```@repl 1
begin
    ActiveVirtualOrbital = new_space(:ActiveVirtualOrbital, "av", "abcdefg");
    ActiveOccupiedOrbital = new_space(:ActiveOccupiedOrbital, "ao", "ijklmno");
    InactiveVirtualOrbital = new_space(:InactiveVirtualOrbital, "iv", "abcdefg");
    InactiveOccupiedOrbital = new_space(:InactiveOccupiedOrbital, "io", "ijklmno");

    set_color(ActiveVirtualOrbital, :cyan)
    set_color(ActiveOccupiedOrbital, :cyan)

    set_color(InactiveVirtualOrbital, :red)
    set_color(InactiveOccupiedOrbital, :red)
end
```

Then we call [`add_space_sum`](@ref) to signal how each space can split into
two smaller spaces

```@repl 1
add_space_sum(ActiveOrbital, InactiveOrbital, GeneralOrbital)
add_space_sum(ActiveVirtualOrbital, InactiveVirtualOrbital, VirtualOrbital)
add_space_sum(ActiveOccupiedOrbital, InactiveOccupiedOrbital, OccupiedOrbital)
add_space_sum(ActiveOccupiedOrbital, ActiveVirtualOrbital, ActiveOrbital)
add_space_sum(InactiveOccupiedOrbital, InactiveVirtualOrbital, InactiveOrbital)
```

Now in addition to this we now have a new relation we want to specify, namely
the following

```
OccupiedOrbital ∩ ActiveOrbital == ActiveOccupiedOrbital
OccupiedOrbital ∩ InactiveOrbital == InactiveOccupiedOrbital
...
```

To achieve this we call the [`add_space_intersection`](@ref) function

```@repl 1
add_space_intersection(ActiveOrbital, OccupiedOrbital, ActiveOccupiedOrbital)
add_space_intersection(ActiveOrbital, VirtualOrbital, ActiveVirtualOrbital)
add_space_intersection(InactiveOrbital, OccupiedOrbital, InactiveOccupiedOrbital)
add_space_intersection(InactiveOrbital, VirtualOrbital, InactiveVirtualOrbital)
```

This is what allows the following simplification to happend

```@repl 1
t = real_tensor("h", 1, 2) * constrain(1 => ActiveOrbital, 2 => InactiveOrbital)
o = E(1, 2) * occupied(1, 2)
t * o
```

where we see that the two indices got constrained to the smaller
`ActiveOccupiedOrbital` and `InactiveOccupiedOrbital` spaces respectively.

We can make some convenience constructor functions

```@repl 1
active(indices...) = constrain(p => ActiveOrbital for p in indices)
inactive(indices...) = constrain(p => InactiveOrbital for p in indices)
aocc(indices...) = constrain(p => ActiveOccupiedOrbital for p in indices)
iocc(indices...) = constrain(p => InactiveOccupiedOrbital for p in indices)
avir(indices...) = constrain(p => ActiveVirtualOrbital for p in indices)
ivir(indices...) = constrain(p => InactiveVirtualOrbital for p in indices)
```

### HF interaction energy

We can now use this to, for example, get an expression for the interaction
energy between the active and inactive spaces at the HF level.

We define the full hamiltonian as

```@repl 1
h = ∑(rsym_tensor("h", 1, 2) * E(1, 2) * electron(1, 2), 1:2)
g = 1//2 * simplify(
    ∑(rsym_tensor("g", 1:4...) * e(1:4...) * electron(1:4...), 1:4)
)
H = h + g
```

which can give us the full HF energy

```@repl 1
E_HF = simplify_heavy(hf_expectation_value(H))
```

We can now also define Hamiltonians acting on only the active and inactive
spaces as

```@repl 1
h_active = ∑(rsym_tensor("h", 1, 2) * E(1, 2) * active(1, 2), 1:2)
g_active = 1//2 * simplify(
    ∑(rsym_tensor("g", 1:4...) * e(1:4...) * active(1:4...), 1:4)
)
h_inactive = ∑(rsym_tensor("h", 1, 2) * E(1, 2) * inactive(1, 2), 1:2)
g_inactive = 1//2 * simplify(
    ∑(rsym_tensor("g", 1:4...) * e(1:4...) * inactive(1:4...), 1:4)
)
H_active = h_active + g_active
H_inactive = h_inactive + g_inactive
```

We can now define the HF energies of the subspaces

```@repl 1
E_active = simplify_heavy(hf_expectation_value(H_active))
E_inactive = simplify_heavy(hf_expectation_value(H_inactive))
```

Then finally we can define an interaction energy by subtracting these off the
full energy

```@repl 1
E_int = simplify_heavy(E_HF - E_active - E_inactive,
    (OccupiedOrbital => (ActiveOccupiedOrbital, InactiveOccupiedOrbital),))
```

Note here the use of the extra argument to the [`simplify_heavy`](@ref)
function. This defaults to splitting general indices into occupied and virtual,
but it was now more useful to split occupied orbitals into active and inactive
parts.


### Multilevel coupled cluster

If we now define the Hamiltonian in terms of the fock matrix as usual for
coupled cluster derivations

```@repl 1
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

Then as an interesting example we can define the T2 cluster operator acting on
active orbitals for one excitation, then for all orbitals for the other

```@repl 1
T2 = 1//2 * ∑(psym_tensor("t", 1:4...) * E(1, 2) * E(3, 4) *
    occupied(2, 4) * virtual(1, 3) * active(1, 2, 3, 4), 1:4) +
    ∑(psym_tensor("t", 1:4...) * E(1, 2) * E(3, 4) *
    occupied(2, 4) * virtual(1, 3) * active(1, 2) * inactive(3, 4), 1:4)
```

We then define the similarity transformed Hamiltonian and its right projection
as usual

```@repl 1
Hbar = simplify(bch(H, T2, 4));
Hbar_ket = simplify(act_on_ket(Hbar, 2));
```

Energy:

```@repl 1
simplify_heavy(act_on_bra(Hbar_ket));
look_for_tensor_replacements(ans,
    make_exchange_transformer("g", "L"));
E0 = ans
```

singles:

```@repl 1
project_biorthogonal(Hbar_ket, E(1, 2));
simplify_heavy(ans, [
    GeneralOrbital => (OccupiedOrbital, VirtualOrbital),
    OccupiedOrbital => (ActiveOccupiedOrbital, InactiveOccupiedOrbital),
    VirtualOrbital => (ActiveVirtualOrbital, InactiveVirtualOrbital),
    GeneralOrbital => (ActiveOrbital, InactiveOrbital),
    ActiveOrbital => (ActiveOccupiedOrbital, ActiveVirtualOrbital),
    InactiveOrbital => (InactiveOccupiedOrbital, InactiveVirtualOrbital),
]);
look_for_tensor_replacements(ans,
    make_exchange_transformer("t", "u"));
look_for_tensor_replacements(ans,
    make_exchange_transformer("g", "L"));
omega_ai = ans
```

Since the different terms here have different constraints for the external
indices, it can be useful to separate this into active/inactive blocks:

```@repl 1
omega_ai_AA = omega_ai * active(1, 2)
omega_ai_AI = omega_ai * active(1) * inactive(2)
omega_ai_IA = omega_ai * active(2) * inactive(1)
omega_ai_II = omega_ai * inactive(1, 2)
```
