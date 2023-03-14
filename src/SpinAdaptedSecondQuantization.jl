module SpinAdaptedSecondQuantization

export SASQ
"""
    SASQ = SpinAdaptedSecondQuantization

Short acronym for accessing internals (like enable_/disable_color)
since the full module name is a bit long.
"""
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

include("hf_expectation_value.jl")
include("ket.jl")
include("wick_theorem.jl")

end # module SpinAdaptedSecondQuantization
