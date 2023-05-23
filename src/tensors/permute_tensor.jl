export permute_tensor

# Permute tensor
# See Eq. 14.4.23 and 14.4.48 in MEST by Helgaker et al.

# Important properties:
#   P_aibj = P_bjai
#   P_aibj F_ai g_bjck = F_ai g_bjck + F_bj g_aick
#   P_aibj F_ai g_bjck = P_aibj F_bj g_aick
#   Has the same symmetries as ParticleSymmetricTensor

struct PermuteTensor <: Tensor
    indices::Vector{Int}
end

function Base.isless(a::PermuteTensor, b::PermuteTensor)
    a.indices < b.indices
end

function Base.isless(a::A, b::B) where {A<:PermuteTensor,B<:Tensor}
    true
end

function Base.isless(a::A, b::B) where {A<:Tensor,B<:PermuteTensor}
    false
end

get_symbol(t::PermuteTensor) = "P"
get_indices(t::PermuteTensor) = t.indices

function sort_permute_indices!(indices)
    sort!(reinterpret(NTuple{2,Int}, indices))
end

function exchange_indices(t::PermuteTensor, mapping)
    new_ind = [exchange_index(i, mapping) for i in t.indices]
    sort_permute_indices!(new_ind)
    PermuteTensor(new_ind)
end

function permute_tensor(indices...)
    ind = collect(indices)
    sort_permute_indices!(ind)
    Expression(PermuteTensor(ind))
end
