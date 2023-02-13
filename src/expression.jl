
struct Expression{T<:Number}
    terms::Vector{Term{T}}

    function Expression(terms::Vector{Term{T}}) where {T<:Number}
        new{T}(sort(terms))
    end
end

function Base.show(io::IO, ex::Expression)
    if isempty(ex.terms)
        throw("Expressions should not be empty,\
but rather include a single zero term")
    end

    print(io, first(ex.terms))

    for t in ex.terms[2:end]
        if t.scalar < 0
            print(io, " - ", new_scalar(t, -t.scalar))
        else
            print(io, " + ", t)
        end
    end
end

function Base.zero(::Type{Expression})
    Expression([zero(Term)])
end

function Expression(s::T) where {T<:Number}
    Expression([Term(
        s,
        MOIndex[],
        KroneckerDelta[],
        Tensor[],
        Operator[]
    )])
end

function Expression(d::KroneckerDelta)
    Expression([Term(
        1,
        MOIndex[],
        [d],
        Tensor[],
        Operator[]
    )])
end

Expression(t::Tensor) = Expression([Term(
    1,
    MOIndex[],
    KroneckerDelta[],
    Tensor[t],
    Operator[]
)])

Expression(o::Operator) = Expression([Term(
    1,
    MOIndex[],
    KroneckerDelta[],
    Tensor[],
    Operator[o]
)])

export summation, ∑

function summation(e::Expression, sum_indices)
    Expression([summation(t, sum_indices) for t in e.terms])
end

# Unicode alias
∑(e, s) = summation(e, s)

function promote_scalar(::Type{T}, e::Expression) where {T<:Number}
    Expression([promote_scalar(T, t) for t in e.terms])
end

function Base.promote(a::Expression{A}, b::Expression{B}) where
{A<:Number,B<:Number}
    (promote_scalar(B, a), promote_scalar(A, b))
end

# Addition:

function Base.:+(a::Expression, b::Expression)
    a, b = promote(a, b)

    Expression([a.terms; b.terms])
end

function Base.:+(a::A, b::Expression) where A <: Number
    Expression(a) + b
end

function Base.:+(a::Expression, b::B) where B <: Number
    a + Expression(b)
end


