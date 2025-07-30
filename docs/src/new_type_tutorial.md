
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
    convenience. See [`Internals`](@ref).

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

Optionally we can implement the functions, but they have reasonable default
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

## Implementing your own operator type

Here as an example we will implement the "hole" excitation operator defined as

``
E^V_{pq} = a_{p\alpha} a^\dagger_{q\alpha} + a_{p\beta} a^\dagger_{q\beta}
=
2 \delta_{pq} - E_{qp}
``

!!! note
    Since the hole excitation operator is expressible in terms of a Kronecker
    delta minus a normal excitation operator, it might be better to just alias
    a constructor for this relation like

    ```@repl 2
    using SpinAdaptedSecondQuantization
    Ev(p, q) = 2 * δ(p, q) - E(q, p)
    Ev(1, 2) * electron(1, 2)
    ```

We start by defining our new type

```@repl 1
struct HoleExcitationOperator <: SASQ.Operator # Must be subtype of Operator
    p::Int
    q::Int
end
```

The set of functions to implement is similar to those needed for implementing
a tensor above.
Required functions are:
- `Base.show` for printing.
- `exchange_indices`
- `get_all_indices`
- `act_on_ket`
- `Base.:(==)`
- `Base.isless`
- `Base.adjoint`

The `Base.show` function takes in a tuple containing the information to
correctly print the names of indices.

```@repl 1
function Base.show(io::IO,
(
    o, constraints, translation
)::Tuple{HoleExcitationOperator,SASQ.Constraints,SASQ.IndexTranslation})
    print(io, "Ev_")
    SASQ.print_mo_index(io, constraints, translation, o.p, o.q)
end
```

```@repl 1
function SASQ.exchange_indices(o::HoleExcitationOperator, mapping)
    HoleExcitationOperator(
        SASQ.exchange_index(o.p, mapping),
        SASQ.exchange_index(o.q, mapping)
    )
end
function SASQ.get_all_indices(o::HoleExcitationOperator)
    (o.p, o.q)
end
function Base.:(==)(a::HoleExcitationOperator, b::HoleExcitationOperator)
    (a.p, a.q) == (b.p, b.q)
end
function Base.isless(a::HoleExcitationOperator, b::HoleExcitationOperator)
    (a.p, a.q) < (b.p, b.q)
end
function Base.adjoint(o::HoleExcitationOperator)
    HoleExcitationOperator(o.q, o.p)
end
function SASQ.act_on_ket(o::HoleExcitationOperator)
    p, q = o.p, o.q

    2 * δ(p, q) * virtual(p, q) -
    E(q, p) * occupied(p) * virtual(q)
end
```

Then there are functions that are implemented for each pair of operator types.
These only needs to be implemented for the pairs you plan on using in the
same expressions. Here we will only be implementing these for the normal
`SingletExcitationOperator`.
- `Base.isless` To define the lexicographic ordering of different operator types
- `reductive_commutator` To define commutation relations

The first function is simple as we only need to decide on an order of the
operator. Say we want to have `E_pq` operators come before `Ev_pq` operators,
then we define the function

```@repl 1
function Base.isless(
    ::Type{HoleExcitationOperator},
    ::Type{SASQ.SingletExcitationOperator})
    false
end
```

next is the `reductive_commutator` which gives the commutation relation between
different operators, giving both the expression for the commutator and whether
it is a commutator or anticommutator. Here both the commutators we want to
implement are commutators, so we return a positive 1 to indicate no sign change
when commuting.

```@repl 1
function SASQ.reductive_commutator(
    a::HoleExcitationOperator,
    b::HoleExcitationOperator)

    p, q = a.p, a.q
    r, s = b.p, b.q

    (1, δ(p, s) * E(q, r) - δ(q, r) * E(s, p))
end
function SASQ.reductive_commutator(
    a::SASQ.SingletExcitationOperator,
    b::HoleExcitationOperator)

    p, q = a.p, a.q
    r, s = b.p, b.q

    (1, δ(p, r) * E(s, q) - δ(q, s) * E(p, r))
end
```

Then finally we make a constructor that produces an `Expression`

```@repl 1
Ev(p, q) = SASQ.Expression(HoleExcitationOperator(p, q))
```

Then we can start using the new operator

```@repl 1
Ev(1, 2) * electron(1, 2)
commutator(Ev(1, 2), E(3, 4)) * electron(1:4...)
∑(fsym_tensor("x", 1, 2, 3, 4) * E(1, 2) * Ev(3, 4) * electron(1:4...), 1:4)
simplify_heavy(act_on_ket(ans))
```
