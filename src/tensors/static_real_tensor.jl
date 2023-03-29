
export static_tensor

# Simplest type of tensor
struct StaticTensor{N} <: Tensor
    symbol::String
    indices::NTuple{N,Int}
end

get_symbol(t::StaticTensor) = t.symbol
get_indices(t::StaticTensor) = t.indices

function exchange_indices(t::StaticTensor, mapping)
    StaticTensor(t.symbol,
        ((exchange_index(i, mapping) for i in t.indices)...,))
end

static_tensor(symbol, indices...) =
    Expression(StaticTensor(symbol, (indices...,)))

