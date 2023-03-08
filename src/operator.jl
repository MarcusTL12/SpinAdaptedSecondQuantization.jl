
abstract type Operator end

include("operators/singlet_excitation_operator.jl")
include("operators/fermion_operator.jl")

function Base.isless(a::A, b::B) where {A<:Operator,B<:Operator}
    !(b < a)
end

function Base.:(==)(::A, ::B) where {A<:Operator,B<:Operator}
    if A == B
        throw("Must implement (==) function for all operator types!")
    else
        false
    end
end

# Implement ordering of new operator types here:

function Base.isless(::FermionOperator, ::SingletExcitationOperator)
    false
end
