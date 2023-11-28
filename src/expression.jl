export bch

struct Expression{T<:Number}
    terms::Vector{Term{T}}

    function Expression(terms::AbstractVector{Term{T}}) where {T<:Number}
        terms = sort(terms)

        # Collect equal termsif isempty(terms)
        if isempty(terms)
            return new{T}([Term(
                zero(T),
                Int[],
                KroneckerDelta[],
                Tensor[],
                Operator[]
            )])
        end

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

function Expression(terms::AbstractVector{Term})
    if isempty(terms)
        return zero(Expression{Int})
    end

    first_term, rest = Iterators.peel(terms)

    T = typeof(first_term.scalar)

    for t in rest
        T = promote_type(T, typeof(t.scalar))
    end

    new_terms = Term{T}[]

    for t in terms
        push!(new_terms, new_scalar(t, T(t.scalar)))
    end

    Expression(new_terms)
end

function Base.show(io::IO, ex::Expression)
    show(io, (ex, IndexTranslation()))
end

function Base.show(io::IO,
    (ex, translation)::Tuple{Expression,IndexTranslation})
    if isempty(ex.terms)
        throw("Expressions should not be empty,\
but rather include a single zero term")
    end

    t, rest = Iterators.peel(ex.terms)

    if t.scalar < 0
        print(io, "- ", (new_scalar(t, -t.scalar), translation))
    else
        print(io, (t, translation))
    end

    for t in rest
        if t.scalar < 0
            print(io, "\n- ", (new_scalar(t, -t.scalar), translation))
        else
            print(io, "\n+ ", (t, translation))
        end
    end
end

function Base.zero(::Type{Expression{T}}) where {T<:Number}
    Expression([zero(Term{T})])
end

function Base.iszero(e::Expression)
    all(iszero, e.terms)
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

export constrain, occupied, virtual

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

function occupied(indices...)
    constrain(p => OccupiedOrbital for p in indices)
end

function virtual(indices...)
    constrain(p => VirtualOrbital for p in indices)
end

function Base.getindex(ex::Expression, i)
    ex.terms[i]
end

function Base.setindex!(ex::Expression, t::Term, i)
    ex.terms[i] = t
end

export scalar_type
function scalar_type(::Expression{T}) where {T<:Number}
    T
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

# Fusion:

function fuse(a::Expression, b::Expression)
    Expression([fuse(t1, t2) for t1 in a.terms for t2 in b.terms])
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

function try_add_constraints(ex::Expression)
    block_begin = 1
    while block_begin <= length(ex.terms)
        block_ident = ex[block_begin]
        block_end = block_begin
        for i in block_begin+1:length(ex.terms)
            if non_constraint_non_scalar_equal(block_ident, ex[i])
                block_end = i
            else
                break
            end
        end

        done = false

        while !done
            done = true

            for i in block_begin:block_end
                for j in i+1:block_end
                    t, did_something = try_add_constraints(ex[i], ex[j])

                    if did_something
                        done = false

                        # n_terms_before = length(ex.terms)

                        if t isa Tuple
                            ex[i] = t[1]
                            ex[j] = t[2]
                        else
                            ex[i] = t
                            deleteat!(ex.terms, j)
                            block_end -= 1
                        end

                        # want to do this, but this is slow
                        # ex = Expression(ex.terms)

                        # So instead do the same work on the subblock
                        sort!(@view ex.terms[block_begin:block_end])
                        k = block_begin
                        while k < block_end
                            while k < block_end &&
                                equal_nonscalar(ex[k], ex[k+1])
                                ex[k] = new_scalar(
                                    ex[k],
                                    ex[k].scalar + ex[k+1].scalar
                                )
                                deleteat!(ex.terms, k + 1)
                                block_end -= 1
                            end
                            k += 1
                        end

                        break
                    end
                end

                if !done
                    break
                end
            end
        end

        block_begin = block_end + 1
    end

    Expression(ex.terms)
end

function simplify_heavy_terms(ex::Expression)
    terms = copy(ex.terms)

    Threads.@threads for i in eachindex(terms)
        terms[i] = simplify_heavy(terms[i])
    end

    Expression(terms)
end

# Unnecessarily expensive simplify for testing on small expressions
export simplify_heavy
function simplify_heavy(ex::Expression)
    done = false
    ex = simplify(ex)
    while !done
        new_ex = try_add_constraints(simplify_heavy_terms(ex))
        done = new_ex == ex
        ex = new_ex
    end
    ex
end

export sort_operators
function sort_operators(ex::Expression)
    Expression([sort_operators(t) for t in ex.terms])
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

# Commutator:
export anticommutator
function anticommutator(a::Expression{A}, b::Expression{B}) where
{A<:Number,B<:Number}
    terms = Term{promote_type(A, B)}[]

    for t1 in a.terms, t2 in b.terms
        append!(terms, anticommutator(t1, t2).terms)
    end

    Expression(terms)
end

function commutator(A::Expression, B::Expression, n::Integer)
    # [A, B]_n
    # [A, B]_3 = [[[A, B], B], B]
    X = A
    for _ = 1:n
        X = commutator(X, B)
    end
    return X
end

function bch(A, B, n)
    # Baker-Campbell-Haussdorf expansion,
    # e^-B A e^B = A + 1/1! [A,B] + 1/2! [[A,B],B] + ... + 1/n! [A,B]_n
    X = A
    Y = A
    for i = 1:n
        Y = 1//i * commutator(Y, B)
        X += Y
    end
    X
end

# Function to express all operators in an expression in terms of
# elementary fermionic/bosinic anihilation and creation operators (if possible)
export convert_to_elementary_operators
function convert_to_elementary_operators(ex::Expression{T}) where {T<:Number}
    terms = Term{T}[]

    for t in ex.terms
        append!(terms, convert_to_elementary_operators(t).terms)
    end

    Expression(terms)
end

function Base.adjoint(ex::Expression)
    Expression([t' for t in ex.terms])
end

function Base.:(^)(ex::Expression{T}, n::Integer) where {T<:Number}
    if n < 0
        throw("Negative powers for expressions are not supported!")
    end

    ex2 = Expression(one(T))
    for _ in 1:n
        ex2 *= ex
    end
    ex2
end
