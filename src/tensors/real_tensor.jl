export real_tensor

# Simplest type of tensor
struct RealTensor <: Tensor
    symbol::String
    indices::Vector{Int}
end

get_symbol(t::RealTensor) = t.symbol
get_indices(t::RealTensor) = t.indices
get_indices_permutations(t::RealTensor) = [t.indices]

function exchange_indices(t::RealTensor, mapping)
    RealTensor(t.symbol, [exchange_index(i, mapping) for i in t.indices])
end

real_tensor(symbol, indices...) =
    Expression(RealTensor(symbol, collect(indices)))
