
# Implementing your own types

In specific cases you might want to implement your own operator and tensor
types. This is a simple matter of implementing the correct set of functions
on the new types. Here is a short tutorial for how to implement your own type.

!!! note
    There is no need to modify the source code of
    `SpinAdaptedSecondQuantization.jl` itself to add new types, as it will
    happily work with whatever types you define in your own input scripts.
    When overloading functions internal to the package, however, the function
    names need to be prepended by `SpinAdaptedSecondQuantization.` to work.
    Since the package name is rather long we export the acronym `SASQ` for
    convencience. See [`Internals`](@ref).

## Implementing your own tensor type

As an example we will implement a tensor where every single index permutation
is equivalent. This might not be very useful, but it is simple to implement.

We start by defining the new tensor type with convenient member variables

```@repl 1
using SpinAdaptedSecondQuantization
struct FullySymmetricTensor <: SASQ.Tensor # Must be subtype to Tensor
    symbol::String
    indices::Vector{Int}
end
```

Next we must implement the following set of functions:
- `get_symbol(t)`, should return the name of the tensor
- `get_indices(t)`, should return ordered list of indices
- `exchange_indices(t, mapping)`, should return tensor after indices are changed
    according to the given mapping
- `reorder_indices(t, permutation)` should return tensor after reordering
    indices according to the given permutation

Optionally we can impement the functions, but they have reasonable default
implementations
- `Base.show(io::IO, t)`, if one wants to print the tensor in a custom way.
- `get_permutations(t)`, needed for [`print_eT_function_generator`](@ref)
    to take advantage of available permutations.

Here we will only implement the required functions.

```@repl 1
SASQ.get_symbol(t::FullySymmetricTensor) = t.symbol
SASQ.get_indices(t::FullySymmetricTensor) = t.indices
function SASQ.exchange_indices(t::FullySymmetricTensor, mapping)
    # Use internal function exchange_index on each index
    new_inds = [SASQ.exchange_index(p, mapping) for p in t.indices]
    sort!(new_inds) # Sort indices as any permutation is equivalent
    FullySymmetricTensor(t.symbol, new_inds)
end
function SASQ.reorder_indices(t::FullySymmetricTensor, permutation)
    new_inds = t.indices[permutation]
    sort!(ind)
    FullySymmetricTensor(t.symbol, new_inds)
end
```

It is also useful to make a constructor that produces an `Expression` containing
that tensor:

```@repl 1
function fsym_tensor(symbol, indices...)
    SASQ.Expression(FullySymmetricTensor(symbol, sort!(collect(indices))))
end
```

Now we can use our new tensor type

```@repl 1
fsym_tensor("x", 2, 3, 2, 7, 5) # Sorts indices
∑(fsym_tensor("x", 1, 2, 3) * E(1, 2) * E(1, 3) * electron(1, 2, 3), 1:3)
act_on_ket(ans)
simplify(ans)
```
