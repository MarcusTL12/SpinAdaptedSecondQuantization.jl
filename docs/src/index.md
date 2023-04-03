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

```@meta
DocTestSetup = quote
  using SpinAdaptedSecondQuantization
  SASQ.disable_color()
end
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
- A dictionary of `orbital constraints`

### Orbital Indices

Most types will contain orbital indices.
For convenience these are represented as integers,
but they are purely symbolic. When printed, they will be given names according
to the semi-standard names for general molecular orbital (MO) indices
`pqrstuvw`. This means that index 1 will be printed as 'p',
2 as 'q' and 8 as 'w'. For indices larger than 8 the names wrap around and
are numbered with subscripts so index 9 will be printed as 'p₁'
and index 78 as 'u₉'.

### Operators

A small variety of operator types is supported, see (REF TO LIST) for full list.

A few examples are:

```@repl 1
E(1, 2)
E(1, 2) * E(3, 4)
e(1, 2, 3, 4)
fermion(1, α)
fermiondag(2, β)
```

### Tensors

The simplest tensor type is the `RealTensor` which has a name and an array of
MO-indices.

Example:

```@repl 1
real_tensor("h", 1, 2)
real_tensor("g", 1, 2, 3, 4)
```

Some other types of tensors are supported with various symmetries in the
indices. See (REF TO LIST) for full list.

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

A term can represent a sum over all values of a specific (or many) MO-indices.

Example:

```@repl 1
a = summation(E(1, 2), [1])
b = ∑(E(1, 2), 1:2) # Unicode
a * b # Automatically renames summation indices to not collide
```

# Orbital Constraints

An important feature is the ability to tell a term whether an MO-index
belongs to a particular subset of all orbitals, specifically whether it belongs
to the occupied or virtual set.

Example:

```@repl 1
occupied(1)
virtual(2)
```

here the constraints are printed explicitly as that is all that exists in
the term, but if the index shows up other places in the term, the indices
are colored to indicate the constraints instead:

```@repl 1
a = E(1, 2) * virtual(1) * occupied(2)
```

Where occupied orbitals are colored green and virtual are colored blue.
If you would rather not have colored prints
(printing to file, no color support in terminal, preferance, ...)
you can easily turn this off:

```@repl 1
SpinAdaptedSecondQuantization.disable_color()
a
SASQ.enable_color() # Can use acronym SASQ instead of full module name
a
```

Constraints will transfer to other indices through Kronecker deltas
and will produce a zero-term if they are unsatisfiable:

```@repl 1
a = δ(1, 2) * occupied(1) # Now q is also occupied
b = δ(2, 3) * virtual(3)
a * b # gives 0 as they cannot be occupied and virtual
```

## Hartree-Fock energy expression

```@repl 1
h = ∑(rsym_tensor("h", 1, 2) * E(1, 2), 1:2)
g = simplify(∑(rsym_tensor("g", 1:4...) * e(1:4...), 1:4))
H = h + g
E_hf = simplify_heavy(act_on_ket(H, 0))
```
