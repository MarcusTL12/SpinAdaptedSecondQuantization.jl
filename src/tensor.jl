
abstract type Tensor end

# These methods should be implemented for all subtypes of Tensor
get_symbol(::T) where {T<:Tensor} =
    throw("get_symbol not implemented for Tensor type $(T)!")
get_indices(::T) where {T<:Tensor} =
    throw("get_indices not implemented for Tensor type $(T)!")
exchange_indices(::T, mapping) where {T<:Tensor} =
    throw("exchange_indices not implemented for Tensor type $(T)!")
get_indices_permutations(::T) where {T<:Tensor} =
    throw("get_indices_permutations not implemented for Tensor type $(T)!")
reorder_indices(::T, permutation) where {T<:Tensor} =
    throw("reorder_indices not implemented for Tensor type $(T)!")

# Base.show is overridable if wanted (typically for cluster amplitudes)
function Base.show(io::IO,
    (t, constraints, translation)::Tuple{Tensor,Constraints,IndexTranslation})
    print(io, get_symbol(t))
    inds = get_indices(t)
    if !isempty(inds)
        print(io, '_')
        for ind in inds
            print_mo_index(io, constraints, translation, ind)
        end
    end
end

function Base.:(==)(a::A, b::B) where {A<:Tensor,B<:Tensor}
    (get_symbol(a), get_indices(a)) == (get_symbol(b), get_indices(b))
end

function Base.isless(a::A, b::B) where {A<:Tensor,B<:Tensor}
    (get_symbol(a), get_indices(a)) < (get_symbol(b), get_indices(b))
end

function get_permutations(t::Tensor)
    [collect(1:length(get_indices(t)))]
end

include("tensors/real_tensor.jl")
include("tensors/psym_tensor.jl")
include("tensors/rsym_tensor.jl")
include("tensors/permute_tensor.jl")
