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

function sort_psym_indices!(indices)
    sort!(reinterpret(NTuple{2,Int}, indices))
end

function exchange_indices(t::ParticleSymmetricTensor, mapping)
    new_ind = [exchange_index(i, mapping) for i in t.indices]
    sort_psym_indices!(new_ind)
    ParticleSymmetricTensor(t.symbol, new_ind)
end

function psym_tensor(symbol, indices...)
    ind = collect(indices)
    @assert iseven(length(ind)) "psym tensors must have an even number of indices"
    sort_psym_indices!(ind)
    Expression(ParticleSymmetricTensor(symbol, ind))
end

function reorder_indices(t::ParticleSymmetricTensor, permutation)
    ind = t.indices[permutation]
    sort_psym_indices!(ind)
    ParticleSymmetricTensor(t.symbol, ind)
end

function get_permutations(t::ParticleSymmetricTensor)
    n = length(get_indices(t)) ÷ 2

    pair_perms = [p.data for p in PermGen(n)]

    function pair_perm_to_perm(pperm)
        perm = Int[]
        for i in pperm
            append!(perm, (2i - 1, 2i))
        end
        perm
    end

    map(pair_perm_to_perm, pair_perms)
end
