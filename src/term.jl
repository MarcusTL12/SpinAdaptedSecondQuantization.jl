
struct Term{T<:Number}
    scalar::T
    deltas::Vector{KroneckerDelta}
    tensors::Vector{Tensor}
    operators::Vector{Operator}
end
