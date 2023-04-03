module SpinAdaptedSecondQuantization

export SASQ
const SASQ = SpinAdaptedSecondQuantization

using DataStructures
using Permutations

"""
    Constraints = SortedDict{Int,Type}

Type alias for container of MO-Index constraints
"""
const Constraints = SortedDict{Int,Type}

include("orbital_spaces.jl")
include("spin.jl")
include("kronecker_delta.jl")
include("operator.jl")
include("tensor.jl")
include("term.jl")
include("expression.jl")

include("ket.jl")
include("wick_theorem.jl")

include("print_code.jl")

end # module SpinAdaptedSecondQuantization
