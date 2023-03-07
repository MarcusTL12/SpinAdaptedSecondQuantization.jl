export rsym_tensor

# Real-Symmetric tensor, chemist notation. (8-fold symmetry)
# (12) = (21)
# (12 34) = (34 12), (12 34) = (21 34), (12 34) = (12 43)

struct RealSymmetricTensor <: Tensor
    symbol::String
    indices::Vector{Int}
end

get_symbol(t::RealSymmetricTensor) = t.symbol
get_indices(t::RealSymmetricTensor) = t.indices

# function get_indices_permutations(t::RealSymmetricTensor)
#     n_pairs = length(t.indices) ÷ 2
#     #n_perm = 2^(n_pairs) * factorial(n_pairs)
#     if n_pairs == 1
#         return [t.indices[[1,2]], t.indices[[2,1]]]
#     elseif n_pairs == 2
#         return [t.indices[[1,2,3,4]], t.indices[[3,4,1,2]],
#                 t.indices[[1,2,4,3]], t.indices[[3,4,2,1]],
#                 t.indices[[2,1,3,4]], t.indices[[4,3,1,2]],
#                 t.indices[[2,1,4,3]], t.indices[[4,3,2,1]]]
#     else
#         throw("not implemented RSymmTensor nPairs = $n_pairs")
#     end
# end

# function sort_rsym_indices(indices)
#     @assert iseven(length(indices))

#     nsets = length(indices) ÷ 2
#     sets = Vector{Vector{MOIndex}}(undef, nsets)
#     for i = 1:nsets
#         sets[i] = [indices[1+2*(i-1)], indices[2+2*(i-1)]]
#     end

#     for i = 1:nsets
#         sort!(sets[i])
#     end
#     sort!(sets)

#     for i = 1:nsets
#         indices[1+2*(i-1)] = sets[i][1]
#         indices[2+2*(i-1)] = sets[i][2]
#     end

#     return indices
# end

function sort_rsym_indices!(indices)
    paired = reinterpret(NTuple{2,Int}, indices)
    for i in eachindex(paired)
        paired[i] = minmax(paired[i]...)
    end
    sort!(paired)
end

function exchange_indices(t::RealSymmetricTensor, mapping)
    new_ind = [exchange_index(i, mapping) for i in t.indices]
    sort_rsym_indices!(new_ind)
    RealSymmetricTensor(t.symbol, new_ind)
end

function rsym_tensor(symbol, indices...)
    ind = collect(indices)
    sort_rsym_indices!(ind)
    Expression(RealSymmetricTensor(symbol, ind))
end
