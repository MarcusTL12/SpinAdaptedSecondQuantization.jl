# Lets you call the constraints as a function instead of explicitly calling get
function (constraints::Constraints)(p::Int)
    get(constraints, p, GeneralOrbital)
end

struct Term{T<:Number}
    scalar::T
    sum_indices::Vector{Int}
    deltas::Vector{KroneckerDelta}
    tensors::Vector{Tensor}
    operators::Vector{Operator}

    # Dict that stores which orbital space each index belongs to
    constraints::Constraints

    function Term(scalar::T, sum_indices, deltas, tensors, operators,
        constraints) where {T<:Number}
        sort!(sum_indices)
        sort!(tensors)

        deltas = compact_deltas(deltas)

        if deltas == 0 || iszero(scalar)
            return new{T}(zero(T), Int[], KroneckerDelta[], Tensor[],
                Operator[], Constraints())
        end

        constraints = copy(constraints)

        # Make mapping from all indices that show up in deltas to the
        # first index in their respective deltas
        nonfirst_delta_indices = Tuple{Int,Int}[]

        # First tighten in the constraints of the first index in each delta
        for d in deltas
            firstind, rest = Iterators.peel(d.indices)
            for p in rest
                s = typeintersect(constraints(firstind), constraints(p))
                if s == Union{}
                    return new{T}(zero(T), Int[], KroneckerDelta[], Tensor[],
                        Operator[], Constraints())
                elseif is_strict_subspace(s, GeneralOrbital)
                    constraints[firstind] = s
                end
                push!(nonfirst_delta_indices, (p, firstind))
            end
        end

        # Then copy those constraints over to all the other indices in the delta
        for (p, q) in nonfirst_delta_indices
            if haskey(constraints, q)
                constraints[p] = constraints[q]
            end
        end

        gen_constraints = Int[]
        for (p, s) in constraints
            if s == GeneralOrbital
                push!(gen_constraints, p)
            end
        end

        for p in gen_constraints
            delete!(constraints, p)
        end

        new{T}(scalar, sum_indices, deltas, tensors, operators, constraints)
    end
end

function Term(scalar::T, sum_indices, deltas, tensors, operators) where
{T<:Number}
    Term(scalar, sum_indices, deltas, tensors, operators, Constraints())
end

Base.copy(t::Term) = Term(
    copy(t.scalar),
    copy(t.sum_indices),
    copy(t.deltas),
    copy(t.tensors),
    copy(t.operators),
    copy(t.constraints)
)

function noop_part(t::Term)
    Term(
        t.scalar,
        t.sum_indices,
        t.deltas,
        t.tensors,
        Operator[],
        t.constraints
    )
end

function Base.zero(::Type{Term{T}}) where {T<:Number}
    Term(zero(T), Int[], KroneckerDelta[], Tensor[], Operator[])
end

function Base.iszero(t::Term)
    iszero(t.scalar)
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
        printscalar(io, t.scalar)
        sep[] = true
    end

    if !isempty(t.sum_indices)
        printsep()
        print(io, "∑_")
        for i in t.sum_indices
            print_mo_index(io, t.constraints, i)
        end
        print(io, '(')
        sep[] = false
    end

    for d in t.deltas
        printsep()
        print(io, t.constraints, d)
    end

    for ten in t.tensors
        printsep()
        print(io, t.constraints, ten)
    end

    for op in t.operators
        printsep()
        print(io, t.constraints, op)
    end

    if !isempty(t.sum_indices)
        print(io, ')')
    end

    constraint_noprint = index_color ? get_non_constraint_indices(t) : Int[]
    constraint_print = [i for (i, _) in t.constraints if i ∉ constraint_noprint]

    if !isempty(constraint_print)
        printsep()

        print(io, "C(")

        isfirst = true

        for i in constraint_print
            if !isfirst
                print(io, ", ")
            end
            print_mo_index(io, i)
            print(io, "∈", getshortname(t.constraints(i)))
            isfirst = false
        end

        print(io, ')')
    end
end

# utility function to "copy" a term but replace the scalar with a new one
function new_scalar(t::Term{T1}, scalar::T2) where {T1<:Number,T2<:Number}
    Term(scalar, t.sum_indices, t.deltas, t.tensors, t.operators, t.constraints)
end

function Base.:-(t::Term)
    new_scalar(t, -t.scalar)
end

function promote_scalar(::Type{T}, t::Term) where {T<:Number}
    new_scalar(t, promote(zero(T), t.scalar)[2])
end

function equal_nonscalar(a::Term, b::Term)
    a.sum_indices == b.sum_indices &&
        a.deltas == b.deltas &&
        a.tensors == b.tensors &&
        a.operators == b.operators &&
        a.constraints == b.constraints
end

function Base.isless(a::Constraints, b::Constraints)
    for (x, y) in zip(a, b)
        if x < y
            return true
        elseif x > y
            return false
        end
    end

    length(a) < length(b)
end

# Exactly how to sort terms is up for debate, but it should be consistent
function Base.isless(a::Term, b::Term)
    (length(a.operators), a.tensors, a.operators, a.deltas, a.sum_indices, a.constraints, b.scalar) <
    (length(b.operators), b.tensors, b.operators, b.deltas, b.sum_indices, b.constraints, a.scalar)
end

function Base.:(==)(a::Term, b::Term)
    (a.tensors, a.operators, a.deltas, a.sum_indices, a.constraints,
        b.scalar) ==
    (b.tensors, b.operators, b.deltas, b.sum_indices, b.constraints,
        a.scalar)
end

function exchange_indices(t::Term{T}, mapping) where {T<:Number}
    if isempty(mapping)
        return t
    end

    t = copy(t)

    for (i, old_ind) in enumerate(t.sum_indices)
        t.sum_indices[i] = exchange_index(old_ind, mapping)
    end

    delete_deltas = Int[]
    for (i, old_delta) in enumerate(t.deltas)
        new_delta = exchange_indices(old_delta, mapping)

        if new_delta isa KroneckerDelta
            t.deltas[i] = new_delta
        elseif new_delta == 1
            push!(delete_deltas, i)
        elseif new_delta == 0
            @warn "Index exchange lead to delta producing zero!"
            return Expression(zero(T))
        end
    end

    for i in reverse!(delete_deltas)
        deleteat!(t.deltas, i)
    end

    for (i, tensor) in enumerate(t.tensors)
        t.tensors[i] = exchange_indices(tensor, mapping)
    end

    for (i, operator) in enumerate(t.operators)
        t.operators[i] = exchange_indices(operator, mapping)
    end

    sort!(t.sum_indices)
    sort!(t.deltas)
    sort!(t.tensors)

    old_constraints = copy(t.constraints)
    empty!(t.constraints)

    for (p, s) in old_constraints
        in_mapping = findfirst(((r, _),) -> r == p, mapping)

        if isnothing(in_mapping)
            t.constraints[p] = s
        else
            t.constraints[mapping[in_mapping][2]] = s
        end
    end

    t
end

function get_non_constraint_indices(t::Term)
    indices = copy(t.sum_indices)

    for d in t.deltas
        append!(indices, d.indices)
    end

    for tensor in t.tensors
        for i in get_indices(tensor)
            push!(indices, i)
        end
    end

    for o in t.operators
        for i in get_all_indices(o)
            push!(indices, i)
        end
    end

    sort!(indices)
    unique!(indices)
end

function get_all_indices(t::Term)
    indices = copy(t.sum_indices)

    for d in t.deltas
        append!(indices, d.indices)
    end

    for tensor in t.tensors
        for i in get_indices(tensor)
            push!(indices, i)
        end
    end

    for o in t.operators
        for i in get_all_indices(o)
            push!(indices, i)
        end
    end

    for (i, _) in t.constraints
        push!(indices, i)
    end

    sort!(indices)
    unique!(indices)
end

# This returns the sum indices of a term
# in the order they show up inside the sum
# The ones that do not show up will come last
function get_sum_indices_ordered(t::Term)
    indices = Int[]

    function add_index(i::Int)
        if i ∉ indices && i ∈ t.sum_indices
            push!(indices, i)
        end
    end

    for d in t.deltas
        for i in d.indices
            add_index(i)
        end
    end

    for o in t.operators
        for i in get_all_indices(o)
            add_index(i)
        end
    end

    for tensor in t.tensors
        for i in get_indices(tensor)
            add_index(i)
        end
    end

    for i in t.sum_indices
        add_index(i)
    end

    indices
end

function enumerate_mo_indices(indices)
    # Examples
    # [p, q, r, s] -> [1, 2, 3, 4]
    # [p, q, r, r] -> [1, 2, 3, 3]
    # [r, r, p, q] -> [1, 1, 2, 3]
    dict = Dict{Int,Int}()
    list = zeros(Int, length(indices))
    counter = 1
    for (n, i) in enumerate(indices)
        if i ∈ keys(dict)
            list[n] = dict[i]
        else
            list[n] = dict[i] = counter
            counter += 1
        end
    end
    return list
end

# These two functions rename summing indices such that there are no
# summing indices that collide with the new indices
function make_space_for_index(t::Term, new_index::Int)
    if new_index ∈ t.sum_indices
        mapping = [new_index => next_free_index(get_all_indices(t))]

        exchange_indices(t, mapping)
    else
        t
    end
end

function make_space_for_indices(t::Term, new_indices)
    indices = get_all_indices(t)
    append!(indices, new_indices)
    sort!(indices)
    unique!(indices)
    mapping = Pair{Int,Int}[]

    for new_index in new_indices
        if new_index ∈ t.sum_indices
            unique_index = next_free_index(indices)
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
        Int[t.sum_indices; sum_indices],
        t.deltas,
        t.tensors,
        t.operators,
        t.constraints
    )
end

# This function reorders the summation indices such that they show up
# in a sorted manner within the sum.
# For example, it will do the conversion:
# ∑_ijab(g_biaj) -> ∑_ijab(g_aibj)
function sort_summation_indices(t::Term)
    if isempty(t.sum_indices)
        return t
    end

    order = get_sum_indices_ordered(t)
    mapping = Pair{Int,Int}[]

    for (pos, ind) in enumerate(sort(order))
        if ind != order[pos]
            push!(mapping, order[pos] => ind)
        end
    end

    exchange_indices(t, mapping)
end

# This function reduces the summation indices to be as small as possible
# when used right before sort_summation_indices it will do the following:
# ∑_ijbd(g_bidj) -> ∑_ijab(g_biaj) -> ∑_ijab(g_aibj)
function lower_summation_indices(t::Term)
    indices = get_all_indices(t)

    mapping = Pair{Int,Int}[]

    for i in reverse(t.sum_indices)
        free_index = next_free_index(indices)

        if free_index < i
            push!(mapping, i => free_index)
            push!(indices, free_index)
        end
    end

    exchange_indices(t, mapping)
end

# This function will look for indices in the tensors and operators
# that show up in Kronecker deltas and exchange them for the
# first index that shows up in that delta.
# Example: δ_pi h_pi E_ip -> δ_pi h_pp E_pp
function lower_delta_indices(t::Term)
    mapping = Pair{Int,Int}[]

    for d in t.deltas
        r, rest = Iterators.peel(d.indices)

        for p in rest
            push!(mapping, p => r)
        end
    end

    new_tensors = [exchange_indices(tensor, mapping) for tensor in t.tensors]
    new_ops = [exchange_indices(o, mapping) for o in t.operators]

    Term(t.scalar, t.sum_indices, t.deltas, new_tensors, new_ops, t.constraints)
end

# This function removes summation indices that show up in kronecker deltas,
# replacing them with the index they would be equal to instead.
# This should be run after `lower_delta_indices`
function simplify_summation_deltas(t::Term)
    t = copy(t)

    done = false

    while !done
        done = true
        for (j, p) in enumerate(t.sum_indices)
            for d in t.deltas
                i = findfirst(==(p), d.indices)
                if !isnothing(i)

                    done = false
                    deleteat!(t.sum_indices, j)

                    if i == 1
                        # If the summation index shows up as the first index
                        # of the delta, then we need to rename all the
                        # occurrences of that index with the next one
                        # in the delta. Example:
                        # ∑_i(δ_ijk h_ip E_iq) -> δ_jk h_jp E_jq

                        t = exchange_indices(t, [p => d.indices[2]])
                    else
                        # Otherwise, it does not show up anywhere else in the
                        # term, so we only need to remove it from the delta

                        t = exchange_indices(t, [p => first(d.indices)])
                    end

                    break
                end
            end

            if done == false
                break
            end
        end
    end

    t
end

function non_constraint_non_scalar_equal(a::Term, b::Term)
    (a.tensors, a.operators, a.deltas, a.sum_indices) ==
    (b.tensors, b.operators, b.deltas, b.sum_indices)
end

# Compares two sets of constraints. If they differ by only one constraint,
# it will return this index. Otherwise, if they are equal, do not contain
# the same set of indices, or differ by more than one constraint, it will
# return nothing
function constraints_equal_but_one(a::Term, b::Term)
    a_indices = get_all_indices(a)
    b_indices = get_all_indices(b)

    if a_indices != b_indices
        return nothing
    end

    p_different = nothing
    for p in a_indices
        s1 = a.constraints(p)
        if b.constraints(p) != s1
            if isnothing(p_different)
                p_different = p
            else
                return nothing
            end
        end
    end

    p_different
end

# TODO: add docs
# TODO: add tests
function try_add_constraints(a::Term, b::Term)
    if !non_constraint_non_scalar_equal(a, b)
        return (a, b), false
    end

    p = constraints_equal_but_one(a, b)
    if isnothing(p)
        return (a, b), false
    end

    s1 = a.constraints(p)
    s2 = b.constraints(p)

    # If we can get one single term by fusing spaces we want to do that
    s12 = add_spaces(s1, s2)
    if a.scalar == b.scalar && !isnothing(s12)
        new_constraints = copy(a.constraints)
        new_constraints[p] = s12
        return Term(a.scalar, a.sum_indices, a.deltas, a.tensors,
            a.operators, new_constraints), true
    end

    if is_strict_subspace(s1, s2)
        s1, s2 = s2, s1
        a, b = b, a
    end

    ds = diff_spaces(s1, s2)
    if !isnothing(ds)
        new_constraints = copy(a.constraints)
        new_constraints[p] = ds

        t1 = Term(a.scalar, a.sum_indices, a.deltas, a.tensors,
            a.operators, new_constraints)

        if iszero(a.scalar + b.scalar)
            return t1, true
        else
            return (t1, new_scalar(b, a.scalar + b.scalar)), true
        end
    end

    (a, b), false
end

function permute_all_sum_indices(t::Term)
    min_t = t
    permuted_inds = copy(t.sum_indices)
    mapping = [p => p for p in permuted_inds]
    for perm in Iterators.drop(PermGen(length(t.sum_indices)), 1)
        copy!(permuted_inds, t.sum_indices)
        permute!(permuted_inds, perm.data)
        for (i, (p, q)) in enumerate(zip(t.sum_indices, permuted_inds))
            mapping[i] = p => q
        end
        min_t = min(min_t, exchange_indices(t, mapping))
    end
    min_t
end

# This function is just a composition of other simplification functions
# in the recomended order to obtain a deterministic simplification of
# the term.
# TODO: This requires extensive testing of determinism and correctness.
# The specific steps and the order of them might need to be adjusted.
function simplify(t::Term)
    t |>
    lower_delta_indices |>
    simplify_summation_deltas |>
    lower_summation_indices |>
    sort_summation_indices
end

# Some operator overloading (Not ment for external use):

function Base.:*(a::A, b::Term{B}) where {A<:Number,B<:Number}
    new_scalar(b, a * b.scalar)
end

# We are assuming commutative scalar multiplication
function Base.:*(a::Term{A}, b::B) where {A<:Number,B<:Number}
    b * a
end

function fuse_constraints!(a::Constraints, b::Constraints)
    for (p, s) in b
        if haskey(a, p)
            if isdisjoint(a[p], s)
                return 0
            else
                a[p] = typeintersect(a[p], s)
            end
        else
            a[p] = s
        end
    end

    1
end

# Multiplying two terms makes sure to rename summation indices such that
# even though they might have some overlap between summation indices, then
# will not be treated as the "same" index.
# Examples:
# ∑_ij(h_ij E_ij) * E_ij = ∑_kl(h_kl E_kl E_ij)
# ∑_i(h_ij) * ∑_i(h_ij) = ∑_ik(h_ij h_kj)
function Base.:*(a::Term{A}, b::Term{B}) where {A<:Number,B<:Number}
    b = make_space_for_indices(b, get_all_indices(a))
    a = make_space_for_indices(a, get_all_indices(b))
    scalar = a.scalar * b.scalar

    sum_indices = [a.sum_indices; b.sum_indices]
    deltas = [a.deltas; b.deltas]
    tensors = Tensor[a.tensors; b.tensors]
    operators = Operator[a.operators; b.operators]

    constraints = copy(a.constraints)
    if fuse_constraints!(constraints, b.constraints) == 0
        return zero(Term{typeof(scalar)})
    end

    Term(scalar, sum_indices, deltas, tensors, operators, constraints)
end

# Like multiplication, but does not make space for summation indices.
function fuse(a::Term, b::Term)
    constraints = copy(a.constraints)
    scalar = a.scalar * b.scalar
    if fuse_constraints!(constraints, b.constraints) == 0
        return zero(Term{typeof(scalar)})
    end

    Term(
        scalar,
        [a.sum_indices; b.sum_indices],
        [a.deltas; b.deltas],
        [a.tensors; b.tensors],
        [a.operators; b.operators],
        constraints
    )
end

# Commutator:

function commutator(a::Term{A}, b::Term{B}) where {A<:Number,B<:Number}
    if isempty(a.operators) || isempty(b.operators)
        return Expression(0)
    end

    b = make_space_for_indices(b, get_all_indices(a))
    a = make_space_for_indices(a, get_all_indices(b))

    terms = Term{promote_type(A, B)}[]

    for i in eachindex(a.operators), j in eachindex(b.operators)
        e = commutator(a.operators[i], b.operators[j])

        lhs = Operator[a.operators[1:i-1]; b.operators[1:j-1]]
        rhs = Operator[b.operators[j+1:end]; a.operators[i+1:end]]

        for t in e.terms
            constraints = copy(a.constraints)
            fuse_constraints!(constraints, t.constraints)
            fuse_constraints!(constraints, b.constraints)

            fused = Term(
                a.scalar * t.scalar * b.scalar,
                Int[a.sum_indices; t.sum_indices; b.sum_indices],
                KroneckerDelta[a.deltas; t.deltas; b.deltas],
                Tensor[a.tensors; t.tensors; b.tensors],
                Operator[lhs; t.operators; rhs],
                constraints
            )

            push!(terms, fused)
        end
    end

    Expression(terms)
end

function commutator_fuse(a::Term{A}, b::Term{B}) where {A<:Number,B<:Number}
    if isempty(a.operators) || isempty(b.operators)
        return Expression(0)
    end

    terms = Term{promote_type(A, B)}[]

    for i in eachindex(a.operators), j in eachindex(b.operators)
        e = commutator(a.operators[i], b.operators[j])

        lhs = Operator[a.operators[1:i-1]; b.operators[1:j-1]]
        rhs = Operator[b.operators[j+1:end]; a.operators[i+1:end]]

        for t in e.terms
            constraints = copy(a.constraints)
            fuse_constraints!(constraints, t.constraints)
            fuse_constraints!(constraints, b.constraints)

            fused = Term(
                a.scalar * t.scalar * b.scalar,
                Int[a.sum_indices; t.sum_indices; b.sum_indices],
                KroneckerDelta[a.deltas; t.deltas; b.deltas],
                Tensor[a.tensors; t.tensors; b.tensors],
                Operator[lhs; t.operators; rhs],
                constraints
            )

            push!(terms, fused)
        end
    end

    Expression(terms)
end