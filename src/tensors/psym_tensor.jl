export psym_tensor

# Particle-Symmetric tensor, chemist notation
# (12 34) = (34 12)
# (12 34 56) = (56 12 34) = (34 56 12) = (56 34 12) = (12 56 34) = (34 12 56)

struct ParticleSymmetricTensor <: Tensor
    symbol::String
    indices::Vector{MOIndex}
end

get_symbol(t::ParticleSymmetricTensor) = t.symbol
get_indices(t::ParticleSymmetricTensor) = t.indices

function sort_psym_indices(indices)
    @assert iseven(length(indices))

    nsets = length(indices) ÷ 2
    sets = Vector{Vector{MOIndex}}(undef, nsets)
    for i = 1:nsets
        sets[i] = [indices[1+2*(i-1)], indices[2+2*(i-1)]]
    end
    sort!(sets)

    for i = 1:nsets
        indices[1+2*(i-1)] = sets[i][1]
        indices[2+2*(i-1)] = sets[i][2]
    end

    return indices
end

function exchange_indices(t::ParticleSymmetricTensor, mapping)
    new_ind = [exchange_index(i, mapping) for i in t.indices]
    sort_ind = sort_psym_indices(new_ind)
    ParticleSymmetricTensor(t.symbol, sort_ind)
end

function psym_tensor(symbol, indices...)
    ind = collect(indices)
    sort_ind = sort_psym_indices(ind)
    Expression(ParticleSymmetricTensor(symbol, sort_ind))
end
