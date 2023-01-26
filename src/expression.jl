
struct Expression{T<:Number}
    terms::Vector{Term{T}}
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
