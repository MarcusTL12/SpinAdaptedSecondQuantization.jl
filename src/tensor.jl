
abstract type Tensor end

# These methods should be implemented for all subtypes of Tensor
get_symbol(::T) where {T<:Tensor} =
    throw("get_symbol not implemented for Tensor type $(T)!")
get_indices(::T) where {T<:Tensor} =
    throw("get_indices not implemented for Tensor type $(T)!")

# Base.show is overridable if wanted (typically for cluster amplitudes)
function Base.show(io::IO, t::T) where {T<:Tensor}
    print(io, get_symbol(t), '_')
    for ind in get_indices(t)
        print(io, ind)
    end
end

function Base.:(==)(a::A, b::B) where {A<:Tensor,B<:Tensor}
    (get_symbol(a), get_indices(a)) == (get_symbol(b), get_indices(b))
end

function Base.isless(a::A, b::B) where {A<:Tensor,B<:Tensor}
    (get_symbol(a), get_indices(a)) < (get_symbol(b), get_indices(b))
end

include("tensors/real_tensor.jl")
