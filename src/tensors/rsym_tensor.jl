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

function reorder_indices(t::RealSymmetricTensor, permutation)
    ind = t.indices[permutation]
    sort_rsym_indices!(ind)
    RealSymmetricTensor(t.symbol, ind)
end

function get_permutations(t::RealSymmetricTensor)
    n = length(get_indices(t)) ÷ 2

    pair_perms = [p.data for p in PermGen(n)]

    function pair_perm_to_perm(pperm, swap_ind)
        perm = Int[]
        for i in pperm
            do_swap = swap_ind & 1 != 0
            swap_ind >>= 1
            append!(perm, do_swap ? (2i, 2i - 1) : (2i - 1, 2i))
        end
        perm
    end

    [pair_perm_to_perm(pp, si) for pp in pair_perms for si in 0:2^n - 1]
end
