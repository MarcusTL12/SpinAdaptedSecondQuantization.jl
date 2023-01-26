
struct Term{T<:Number}
    scalar::T
    sum_indices::Vector{MOIndex}
    deltas::Vector{KroneckerDelta}
    tensors::Vector{Tensor}
    operators::Vector{Operator}

    function Term(scalar::T, sum_indices, deltas, tensors, operators) where
    {T<:Number}
        sort!(sum_indices)
        sort!(deltas)
        sort!(tensors)

        new{T}(scalar, sum_indices, deltas, tensors, operators)
    end
end

Base.copy(t::Term) = Term(
    copy(t.scalar),
    copy(t.sum_indices),
    copy(t.deltas),
    copy(t.tensors),
    copy(t.operators)
)

function Base.zero(::Type{Term})
    Term(0, MOIndex[], KroneckerDelta[], Tensor[], Operator[])
end

function printscalar(io::IO, s::T) where {T<:Number}
    print(io, s)
end

function printscalar(io::IO, s::Rational{T}) where {T}
    if isone(denominator(s))
        print(io, numerator(s))
    else
        print(io, numerator(s), "/", denominator(s))
    end
end

function Base.show(io::IO, t::Term{T}) where {T<:Number}
    sep = Ref(false)

    function printsep()
        if sep[]
            print(io, ' ')
        end
        sep[] = true
    end

    all_nonscalar_empty = isempty(t.sum_indices) && isempty(t.deltas) &&
                          isempty(t.tensors) && isempty(t.operators)

    if !isone(t.scalar)
        if isone(-t.scalar)
            print(io, '-')
        else
            printsep()
            printscalar(io, t.scalar)
        end
    elseif all_nonscalar_empty
        print(io, t.scalar)
    end

    if !isempty(t.sum_indices)
        printsep()
        print(io, "∑_")
        for i in t.sum_indices
            print(io, i)
        end
        print(io, '(')
        sep[] = false
    end

    for d in t.deltas
        printsep()
        print(io, d)
    end

    for ten in t.tensors
        printsep()
        print(io, ten)
    end

    for op in t.operators
        printsep()
        print(io, op)
    end

    if !isempty(t.sum_indices)
        print(io, ')')
    end
end

# utility function to "copy" a term but replace the scalar with a new one
function new_scalar(t::Term{T1}, scalar::T2) where {T1<:Number,T2<:Number}
    Term(scalar, t.sum_indices, t.deltas, t.tensors, t.operators)
end

function exchange_indices(t::Term{T}, mapping) where {T<:Number}
    if isempty(mapping)
        return t
    end

    t = copy(t)

    for (i, old_ind) in enumerate(t.sum_indices)
        t.sum_indices[i] = exchange_index(old_ind, mapping)
    end

    for (i, old_delta) in enumerate(t.deltas)
        new_delta = KroneckerDelta(
            exchange_index(old_delta.p, mapping),
            exchange_index(old_delta.q, mapping)
        )

        if new_delta isa KroneckerDelta
            t.deltas[i] = new_delta
        else
            @warn "Index exchange lead to delta producing zero!"
            return Expression(zero(T))
        end
    end

    for (i, tensor) in enumerate(t.tensors)
        t.tensors[i] = exchange_indices(tensor, mapping)
    end

    for (i, operator) in enumerate(t.operators)
        t.operators[i] = exchange_indices(operator, mapping)
    end

    t
end

function get_all_indices(t::Term)
    indices = copy(t.sum_indices)

    function add_index(i::MOIndex)
        if i ∉ indices
            push!(indices, i)
        end
    end

    for d in t.deltas
        add_index(d.p)
        add_index(d.q)
    end

    for tensor in t.tensors
        for i in get_indices(tensor)
            add_index(i)
        end
    end

    for o in t.operators
        for i in get_all_indices(o)
            add_index(i)
        end
    end

    sort!(indices)
end

# These two functions rename summing indices such that there are no
# summing indices that collide with the new indices
function make_space_for_indices(t::Term, new_index::MOIndex)
    if new_index ∈ t.sum_indices
        mapping = [new_index => next_free_index(get_all_indices(t), new_index)]

        exchange_indices(t, mapping)
    else
        t
    end
end

function make_space_for_indices(t::Term, new_indices)
    indices = get_all_indices(t)
    mapping = Pair{MOIndex,MOIndex}[]

    for new_index in new_indices
        if new_index ∈ t.sum_indices
            unique_index = next_free_index(indices, new_index)
            push!(indices, unique_index)
            push!(mapping, new_index => unique_index)
        end
    end

    exchange_indices(t, mapping)
end

function summation(t::Term, sum_indices)
    t = make_space_for_indices(t, sum_indices)
    Term(
        t.scalar,
        MOIndex[t.sum_indices; sum_indices],
        t.deltas,
        t.tensors,
        t.operators
    )
end
