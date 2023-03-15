
"""
    Operator

Abstact operator type which all concrete operator types (i.e. E_pq, a†_p)
must extend.
"""
abstract type Operator end

include("operators/singlet_excitation_operator.jl")
include("operators/fermion_operator.jl")

include("operators/sorting.jl")
include("operators/commutation_relations.jl")

# Required methods for all Operators to overload:

"""
    Base.print(::IO, ::Constraints, ::Operator)

All typed extending `Operator` must overload this function.
"""
function Base.print(::IO, ::Constraints, ::Operator)
    throw("All operator types must overload Base.print(io, constraints, o)")
end

"""
    exchange_indices(::Operator, mapping)

All typed extending `Operator` must overload this function.
It should return a new operator where all indices have been exchanged
according to the `mapping`.
"""
function exchange_indices(::Operator, mapping)
    throw("All operator types must overload exchange_indices")
end

"""
    get_all_indices(::Operator)

All typed extending `Operator` must overload this function.
It should return an iteratable over all the MO-indices contained in the operator
in order.
"""
function get_all_indices(::Operator)
    throw("All operator types must overload get_all_indices")
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
should be added to the src/operators/sorting.jl file. These should be dummy
routines that simply return true or false, and do not depend on the contents
of the operators.

!!! note
    Only one order needs to be implemented, meaning if `isless(::A, ::B)`
    is implemented, there is no need to also implement `isless(::B, ::A)`
"""
function Base.isless(a::A, b::B) where {A<:Operator,B<:Operator}
    !(b < a)
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

Defines how operators commute; TODO: complete this
"""
function reductive_commutator(a::Operator, b::Operator)
    Γ, c = reductive_commutator(b, a)

    (Γ, -Γ * c)
end
