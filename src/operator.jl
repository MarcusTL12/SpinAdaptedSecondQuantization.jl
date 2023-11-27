
"""
    Operator

Abstact operator type which all concrete operator types (i.e. E_pq, a†_p)
must extend.
"""
abstract type Operator end

include("operators/singlet_excitation_operator.jl")
include("operators/triplet_excitation_operator.jl")
include("operators/fermion_operator.jl")
include("operators/boson_operator.jl")

include("operators/sorting.jl")
include("operators/commutation_relations.jl")

# Required methods for all Operators to overload:

"""
    Base.show(::IO, ::Tuple{Operator,Constraints,IndexTranslation})

All typed extending `Operator` must overload this function.
"""
function Base.show(::IO, ::Tuple{Operator,Constraints,IndexTranslation})
    throw("Base.show(::IO, ::Tuple{$O,Constraints,IndexTranslation}) \
    not implemented!")
end

"""
    exchange_indices(::Operator, mapping)

All typed extending `Operator` must overload this function.
It should return a new operator where all indices have been exchanged
according to the `mapping`.
"""
function exchange_indices(::Operator, mapping)
    throw("exchange_indices(::$O, mapping) not implemented!")
end

"""
    get_all_indices(::Operator)

All typed extending `Operator` must overload this function.
It should return an iteratable over all the MO-indices contained in the operator
in order.
"""
function get_all_indices(::Operator)
    throw("get_all_indices(::$O) not implemented!")
end

"""
    act_on_ket(::Operator)

All typed extending `Operator` must overload this function.
It should return the expression you end up with after acting the operator on
a normal RHF ket.

Example:

``E_{p q} |HF⟩ = (2 δ_{p q} C(p∈O,q∈O) + E_{p q} C(p∈V,q∈O)) |HF⟩``
"""
function act_on_ket(::Operator)
    throw("act_on_ket(::$O) not implemented!")
end

# Generic fallback functions to reduce number of required functions to overload

"""
    Base.isless(::Operator, ::Operator)

All operator types must implement how they should be sorted among themselves
(`isless(::O, ::O)`) and with respect to all other operator types that it is
expected to be used together with (`isless(::A, ::B)`).

The implementation for how it sorts among the same operator type should be
implemented in the definition file for that operator type.

The implementation for how the type sorts compared to other operator types
should be added to the `src/operators/sorting.jl` file. These should be dummy
routines that simply return true or false, and do not depend on the contents
of the operators.

!!! note
    Only one order needs to be implemented, meaning if `isless(::A, ::B)`
    is implemented, there is no need to also implement `isless(::B, ::A)`
"""
function Base.isless(::Type{A}, ::Type{B}) where {A<:Operator,B<:Operator}
    !(B < A)
end

"""
    Base.:==(::Operator, ::Operator)

This has to be explicitly implemented for each operator type.
"""
function Base.:(==)(::A, ::B) where {A<:Operator,B<:Operator}
    if A == B
        throw("Must implement (==) function for all operator types!")
    else
        false
    end
end

"""
    reductive_commutator(::Operator, ::Operator)

Defines how operators commute. This will act as either a commutator or
anticommutator depending on which reduces the rank of the resulting operators.

Returns a `Tuple{Int,Expression}` which gives the type of commutator used
and the resulting expression. The type is signaled by returning 1 or -1 as
the first element of the tuple. This represents the sign change by commutation
which is 1 for a normal commutator ``([A, B] = 0 ⇒ AB == BA)`` and -1 for an
anticommutator ``([A, B]_+ = 0 ⇒ AB = -BA)``

This should be implemented for all pairs of operators types that are expected
to show up in the same expression. The implementations should be put in the
`src/operators/commutation_relations.jl` file.

!!! note
    Only one order needs to be implemented, meaning if
    `reductive_commutator(::A, ::B)` is implemented,
    there is no need to also implement `reductive_commutator(::B, ::A)`
"""
function reductive_commutator(a::Operator, b::Operator)
    Γ, c = reductive_commutator(b, a)

    (Γ, -Γ * c)
end
