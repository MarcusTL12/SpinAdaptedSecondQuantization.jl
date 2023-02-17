
abstract type Operator end

function Base.isless(a::O, b::O) where O <: Operator
    get_all_indices(a) < get_all_indices(b)
end

# Implement ordering of new operator types here:

include("operators/singlet_excitation_operator.jl")
include("operators/fermion_operator.jl")