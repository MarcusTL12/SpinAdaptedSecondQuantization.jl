
abstract type Operator end

include("operators/singlet_excitation_operator.jl")
include("operators/fermion_operator.jl")

include("operators/sorting.jl")
include("operators/commutation_relations.jl")

# Generic fallback functions to reduce number of required functions to overload

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

# Swap so only one function has to be implemented per pair of Operator types
# [A, B]_± = ±[B, A]
function reductive_commutator(a::A, b::B) where {A<:Operator,B<:Operator}
    Γ, c = reductive_commutator(b, a)

    (Γ, -Γ * c)
end
