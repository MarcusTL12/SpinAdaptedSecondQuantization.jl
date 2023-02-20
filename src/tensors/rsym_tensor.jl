export rsym_tensor

# Real-Symmetric tensor, chemist notation. (8-fold symmetry)
# (12) = (21)
# (12 34) = (34 12), (12 34) = (21 34), (12 34) = (12 43)

struct RealSymmetricTensor <: Tensor
    symbol::String
    indices::Vector{MOIndex}
end

get_symbol(t::RealSymmetricTensor) = t.symbol
get_indices(t::RealSymmetricTensor) = t.indices

function sort_rsym_indices(indices)
    @assert iseven(length(indices))

    nsets = length(indices) ÷ 2
    sets = Vector{Vector{MOIndex}}(undef, nsets)
    for i = 1:nsets
        sets[i] = [indices[1+2*(i-1)], indices[2+2*(i-1)]]
    end

    for i = 1:nsets
        sort!(sets[i])
    end
    sort!(sets)

    for i = 1:nsets
        indices[1+2*(i-1)] = sets[i][1]
        indices[2+2*(i-1)] = sets[i][2]
    end

    return indices
end

function exchange_indices(t::RealSymmetricTensor, mapping)
    new_ind = [exchange_index(i, mapping) for i in t.indices]
    sort_ind = sort_rsym_indices(new_ind)
    ParticleSymmetricTensor(t.symbol, sort_ind)
end

function rsym_tensor(symbol, indices...)
    ind = collect(indices)
    sort_ind = sort_rsym_indices(ind)
    Expression(RealSymmetricTensor(symbol, sort_ind))
end
