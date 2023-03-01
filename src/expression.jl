
struct Expression{T<:Number}
    terms::Vector{Term{T}}

    function Expression(terms::Vector{Term{T}}) where {T<:Number}
        terms = sort(terms)

        # Collect equal terms

        first_term, rest = Iterators.peel(terms)
        terms = Term{T}[first_term]

        for t in rest
            if equal_nonscalar(last(terms), t)
                terms[end] = new_scalar(
                    terms[end], terms[end].scalar + t.scalar
                )
            else
                push!(terms, t)
            end
        end

        filter!(!iszero, terms)

        if isempty(terms)
            new{T}([Term(
                zero(T),
                Int[],
                KroneckerDelta[],
                Tensor[],
                Operator[]
            )])
        else
            new{T}(terms)
        end
    end
end

function Base.show(io::IO, ex::Expression)
    if isempty(ex.terms)
        throw("Expressions should not be empty,\
but rather include a single zero term")
    end

    t, rest = Iterators.peel(ex.terms)

    if t.scalar < 0
        println(io, "- ", new_scalar(t, -t.scalar))
    else
        println(io, t)
    end

    for t in rest
        if t.scalar < 0
            println(io, "- ", new_scalar(t, -t.scalar))
        else
            println(io, "+ ", t)
        end
    end
end

function Base.zero(::Type{Expression{T}}) where {T<:Number}
    Expression([zero(Term{T})])
end

function Expression(s::T) where {T<:Number}
    Expression([Term(
        s,
        Int[],
        KroneckerDelta[],
        Tensor[],
        Operator[]
    )])
end

function Expression(d::KroneckerDelta)
    Expression([Term(
        1,
        Int[],
        [d],
        Tensor[],
        Operator[]
    )])
end

Expression(t::Tensor) = Expression([Term(
    1,
    Int[],
    KroneckerDelta[],
    Tensor[t],
    Operator[]
)])

Expression(o::Operator) = Expression([Term(
    1,
    Int[],
    KroneckerDelta[],
    Tensor[],
    Operator[o]
)])

Expression(ops::Vector{Operator}) = Expression([Term(
    1,
    Int[],
    KroneckerDelta[],
    Tensor[],
    ops
)])

export constrain

function constrain(constraints...)
    Expression([Term(
        1,
        Int[],
        KroneckerDelta[],
        Tensor[],
        Operator[],
        Constraints(constraints...)
    )])
end

function Base.getindex(ex::Expression, i)
    ex.terms[i]
end

function Base.setindex!(ex::Expression, t::Term, i)
    ex.terms[i] = t
end

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

function Base.:+(a::A, b::Expression) where {A<:Number}
    Expression(a) + b
end

function Base.:+(a::Expression, b::B) where {B<:Number}
    a + Expression(b)
end

# Subtraction:

function Base.:-(e::Expression)
    Expression([-t for t in e.terms])
end

function Base.:-(a::Expression, b)
    a + -b
end

# Multiplication:

function Base.:*(a::Expression, b::B) where {B<:Number}
    Expression([t * b for t in a.terms])
end

function Base.:*(a::A, b::Expression) where {A<:Number}
    b * a
end

function Base.:/(a::Expression, b::B) where {B<:Number}
    a * inv(b)
end

function Base.://(a::Expression, b::B) where {B<:Number}
    a * (1 // b)
end

function Base.:*(a::Expression, b::Expression)
    Expression([t1 * t2 for t1 in a.terms for t2 in b.terms])
end

function Base.:(==)(a::Expression, b::Expression)
    a.terms == b.terms
end

# Simplification:

export simplify
function simplify(ex::Expression)
    ex |>
    simplify_terms |>
    try_add_constraints |>
    simplify_terms
end

function simplify_terms(ex::Expression)
    Expression([simplify(t) for t in ex.terms])
end

# TODO: group terms into equal significant parts
function try_add_constraints(ex::Expression)
    done = false

    while !done
        done = true

        for i in eachindex(ex.terms)
            for j in i+1:length(ex.terms)
                if !non_constraint_non_scalar_equal(ex[i], ex[j])
                    break
                end

                t, did_something = try_add_constraints(ex[i], ex[j])

                if did_something
                    done = false

                    if t isa Tuple
                        ex[i] = t[1]
                        ex[j] = t[2]
                    else
                        ex[i] = t
                        deleteat!(ex.terms, j)
                    end

                    ex = Expression(ex.terms)

                    break
                end
            end

            if !done
                break
            end
        end
    end

    ex
end

# Commutator:

export commutator
function commutator(a::Expression{A}, b::Expression{B}) where
{A<:Number,B<:Number}
    terms = Term{promote_type(A, B)}[]

    for t1 in a.terms, t2 in b.terms
        append!(terms, commutator(t1, t2).terms)
    end

    Expression(terms)
end
