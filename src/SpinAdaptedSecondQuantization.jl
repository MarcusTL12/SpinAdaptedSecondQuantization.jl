module SpinAdaptedSecondQuantization

using DataStructures

const Constraints = SortedDict{Int,Type}
include("orbital_spaces.jl")
include("kronecker_delta.jl")
include("operator.jl")
include("tensor.jl")
include("term.jl")
include("expression.jl")

include("hf_expectation_value.jl")
include("ket.jl")

end # module SpinAdaptedSecondQuantization
