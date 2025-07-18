# SpinAdaptedSecondQuantization.jl

[SpinAdaptedSecondQuantization.jl]
(https://github.com/MarcusTL12/SpinAdaptedSecondQuantization.jl)
is a Julia package for doing symbolic second quantization mainly targeted at
quantum chemistry methods, such as post Hartree-Fock methods.

To install the package, run

```julia
(v1.8) pkg> add SpinAdaptedSecondQuantization
```

Julia version 1.8 or higher is required.

The package can then be loaded

```julia
using SpinAdaptedSecondQuantization
```

```@setup 1
using SpinAdaptedSecondQuantization
```

## Basic Structure

The core type of the package is the `Expression` type. This contains an array of
`Term` types. A term is a product that can contain a combination of various
elements. These are:

- A `scalar` multiplier
- An array of `operators`
- An array of `tensors`
- An array of `Kronecker deltas`
- An array of `summation indices`
- A dictionary of `index constraints`

### Indices

Most types will contain indices.
For convenience these are represented as integers,
but they are purely symbolic. When constrained to a certain index space
they will be printed with names according
to the semi-standard names for general molecular orbital (MO) indices
`pqrstuv`, `abcdefg` and `ijklmno`.

By default indices are unconstrained and will be printed as subscript indices
`₁₂...`

```@repl 1
E(1, 2)
```

Usually, however, one wants to constrain indices to a certain subspace of
indices. By default the available spaces are `GeneralOrbital`, `VirtualOrbital`
and `OccupiedOrbital` which can be used for example like:

```@repl 1
E(1, 2) * constrain(1 => VirtualOrbital, 2 => OccupiedOrbital)
```

The shorthand functions `electron(...)`, `virtual(...)` and `occupied(...)`
can also be used

```@repl 1
E(1, 2) * electron(1, 2)
E(1, 2) * E(3, 4) * virtual(1, 3) * occupied(2, 4)
```

### Operators

A small variety of operator types is supplied by defualt listed in the following
table.

| Operator Type             | constructor(s)                  |
| ------------------------- | ------------------------------- |
| SingletExcitationOperator | E(p, q), e(p, q, r, s)          |
| FermionOperator           | fermion(p, σ), fermiondag(p, σ) |
| BosonOperator             | boson(), bosondag()             |
| TripletExcitationOperator | τ(p, q)                         |

A few examples are:

```@repl 1
E(1, 2) * electron(1, 2)
E(1, 2) * E(3, 4) * virtual(1, 3) * occupied(2, 4)
e(1, 2, 3, 4) * electron(1, 2, 3, 4)
fermion(1, α) * occupied(1)
fermiondag(2, β) * virtual(2)
bosondag() * boson()
τ(1, 2) * electron(1, 2)
```

### Tensors

The simplest tensor type is the `RealTensor` which has a name and an array of
indices.

Example:

```@repl 1
real_tensor("h", 1, 2) * electron(1, 2)
real_tensor("g", 1, 2, 3, 4) * electron(1, 2, 3, 4)
```

Some other types of tensors are supported with various symmetries in the
indices. Supported tensor types are listed below

| Tensor type             | constructor               | symmetries                        |
| ----------------------- | ------------------------- | --------------------------------- |
| RealTensor              | real_tensor(name, p, ...) | No symmetries                     |
| ParticleSymmetricTensor | psym_tensor(name, p, ...) | g₁₂₃₄ <-> g₃₄₁₂                   |
| RealSymmetricTensor     | rsym_tensor(name, p, ...) | g₁₂₃₄ <-> g₂₁₃₄ <-> g₄₃₂₁ <-> ... |

`rsym_tensor` i typically used when assuming full real orbtial symmetry of
integrals such as 2-fold symmetry of h_pq and 8-fold symmetry of g_pqrs.
the `psym_tensor` is useful when doing coupled cluster theory, where symmetries
within index pairs does not exist because the integrals are T1 transformed.

Examples of symmetric tensors:

```@repl 1
rsym_tensor("g", 4, 3, 1, 2) * electron(1, 2, 3, 4)
psym_tensor("g", 4, 3, 1, 2) * electron(1, 2, 3, 4)
```

### Kronecker Deltas

Kronecker deltas constrain indices to be equal or else the term would be zero.
They can have two or more indices each.

Example:

```@repl 1
delta(1, 2)
δ(1, 2) # Unicode version equivalent to the one above
d1 = δ(1, 2)
d2 = δ(2, 3)
d1 * d2 # Compacts delta expression since all are equal
```

### Summation Indices

A term can represent a sum over all values of a specific (or many) indices.
A sum can be constructed using the `summation` function, or the unicode aliases
`∑` or `Σ`.

Example:

```@repl 1
a = summation(E(1, 2) * electron(1, 2), [1])
b = ∑(E(1, 2) * electron(1, 2), 1:2) # Unicode
a * b # Automatically renames summation indices to not collide
```

## Hartree-Fock energy expression

As a simple example of use of the package, here is how to derive the HF energy
expression.

```@repl 1
h = ∑(real_tensor("h", 1, 2) * E(1, 2) * electron(1, 2), 1:2)
g = 1//2 * simplify(
    ∑(psym_tensor("g", 1:4...) * e(1:4...) * electron(1:4...), 1:4)
)
H = h + g + real_tensor("h_nuc")
E_hf = simplify_heavy(hf_expectation_value(H))
```
