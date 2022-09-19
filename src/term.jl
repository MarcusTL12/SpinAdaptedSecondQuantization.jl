
struct Term{T<:Number}
    scalar::T
    deltas::Vector{KroeneckerDelta}
    tensors::Vector{Tensor}
    operators::Vector{Operator}
end
