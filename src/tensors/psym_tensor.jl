export psym_tensor

# Particle-Symmetric tensor, chemist notation
# (12 34) = (34 12)
# (12 34 56) = (56 12 34) = (34 56 12) = (56 34 12) = (12 56 34) = (34 12 56)

struct ParticleSymmetricTensor <: Tensor
    symbol::String
    indices::Vector{Int}
end

get_symbol(t::ParticleSymmetricTensor) = t.symbol
get_indices(t::ParticleSymmetricTensor) = t.indices

# function get_indices_permutations(t::ParticleSymmetricTensor)
#     n_pairs = length(t.indices) ÷ 2
#     if n_pairs == 1
#         return [t.indices]
#     elseif n_pairs == 2
#         ind = copy(t.indices)
#         return [t.indices, permute!(ind, [3,4,1,2])]
#     else
#         throw("not implemented PSymmTensor nPairs = $n_pairs")
#     end
# end

function sort_psym_indices!(indices)
    @assert iseven(length(indices))

    indices_paired = reinterpret(NTuple{2,Int}, indices)
    sort!(indices_paired)
end

function exchange_indices(t::ParticleSymmetricTensor, mapping)
    new_ind = [exchange_index(i, mapping) for i in t.indices]
    sort_psym_indices!(new_ind)
    ParticleSymmetricTensor(t.symbol, new_ind)
end

function psym_tensor(symbol, indices...)
    ind = collect(indices)
    sort_psym_indices!(ind)
    Expression(ParticleSymmetricTensor(symbol, ind))
end
