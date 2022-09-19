
# Simplest type of tensor
struct RealTensor <: Tensor
    name::String
    indices::Vector{MOIndex}
end
