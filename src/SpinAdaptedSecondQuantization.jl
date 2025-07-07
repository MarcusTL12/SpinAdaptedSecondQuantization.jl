module SpinAdaptedSecondQuantization

export SASQ
const SASQ = SpinAdaptedSecondQuantization

using DataStructures
using Permutations

include("index_spaces.jl")
include("index_space_definitions.jl")
include("spin.jl")
include("kronecker_delta.jl")
include("operator.jl")
include("tensor.jl")
include("term.jl")
include("expression.jl")

include("ket.jl")
include("wick_theorem.jl")

include("biorthonormal.jl")

include("code_generation.jl")

include("tensor_replacements.jl")

include("desymmetrization.jl")

include("precompile.jl")

end # module SpinAdaptedSecondQuantization
