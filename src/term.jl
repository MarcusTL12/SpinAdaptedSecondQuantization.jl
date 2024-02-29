struct Term{T<:Number}
    scalar::T
    sum_indices::Vector{Int}
    deltas::Vector{KroneckerDelta}
    tensors::Vector{Tensor}
    operators::Vector{Operator}

    # Dict that stores which orbital space each index belongs to
    constraints::Constraints

    # Tag to say that it is unnessecary to try to further simplify this term
    max_simplified::Bool

    function Term(scalar::T, sum_indices, deltas, tensors, operators,
        constraints, max_simplified, nocheck) where {T<:Number}
        if nocheck
            new{T}(scalar, sum_indices, deltas, tensors, operators,
                constraints, max_simplified)
        else
            Term(scalar, sum_indices, deltas, tensors, operators,
                constraints, max_simplified)
        end
    end

    function Term(scalar::T, sum_indices, deltas, tensors, operators,
        constraints, max_simplified) where {T<:Number}
        sort!(sum_indices)
        sort!(tensors)

        deltas = compact_deltas(deltas)

        if deltas == 0 || iszero(scalar)
            return new{T}(zero(T), Int[], KroneckerDelta[], Tensor[],
                Operator[], Constraints(), true)
        end

        constraints = copy(constraints)

        # Make mapping from all indices that show up in deltas to the
        # first index in their respective deltas
        nonfirst_delta_indices = Tuple{Int,Int}[]

        # First tighten in the constraints of the first index in each delta
        for d in deltas
            firstind, rest = Iterators.peel(d.indices)
            for p in rest
                s = intersect(constraints(firstind), constraints(p))
                if isnothing(s)
                    return new{T}(zero(T), Int[], KroneckerDelta[], Tensor[],
                        Operator[], Constraints(), true)
                elseif s != GeneralIndex
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

        new{T}(
            scalar,
            sum_indices,
            deltas,
            tensors,
            operators,
            constraints,
            max_simplified
        )
    end
end

function Term(scalar::T, sum_indices, deltas, tensors, operators,
    constraints) where {T<:Number}
    Term(scalar, sum_indices, deltas, tensors, operators, constraints, false)
end

function Term(scalar::T, sum_indices, deltas, tensors, operators) where
{T<:Number}
    Term(scalar, sum_indices, deltas, tensors, operators, Constraints())
end

Base.copy(t::Term, max_simplified=t.max_simplified) = Term(
    copy(t.scalar),
    copy(t.sum_indices),
    copy(t.deltas),
    copy(t.tensors),
    copy(t.operators),
    copy(t.constraints),
    max_simplified,
    true
)

function noop_part(t::Term)
    Term(
        t.scalar,
        t.sum_indices,
        t.deltas,
        t.tensors,
        Operator[],
        t.constraints,
        t.max_simplified,
        true
    )
end

function new_constraints(t::Term, constraints::Constraints)
    Term(
        t.scalar,
        t.sum_indices,
        t.deltas,
        t.tensors,
        t.operators,
        constraints,
        false,
        true
    )
end

function Base.zero(::Type{Term{T}}) where {T<:Number}
    Term(zero(T), Int[], KroneckerDelta[], Tensor[], Operator[],
        Constraints(), true, true)
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

# TODO: Make this more subspace agnostic
function update_index_translation(t::Term, translation::IndexTranslation)
    translation = copy(translation)

    ex_inds = get_external_indices(t)

    function find_first_free!(s)
        for i in Iterators.countfrom(1)
            if i ∉ s
                push!(s, i)
                return i
            end
        end
    end

    if do_index_translation
        seen_inds = Dict{Int,Set{Int}}()

        for (p, (S, q)) in translation
            if p ∈ ex_inds
                push!(get!(seen_inds, S, Set()), q)
            end
        end

        for p in ex_inds
            if !haskey(translation, p)
                S = t.constraints(p)
                translation[p] = (S, find_first_free!(get!(seen_inds, S, Set())))
            end
        end

        for p in t.sum_indices
            S = t.constraints(p)
            translation[p] = (S, find_first_free!(get!(seen_inds, S, Set())))
        end
    end

    translation
end

function Base.show(io::IO, (t, translation)::Tuple{Term,IndexTranslation})
    sep = Ref(false)

    function printsep()
        if sep[]
            print(io, ' ')
        end
        sep[] = true
    end

    ex_inds = get_external_indices(t)

    for (p, (S, _)) in translation
        Sc = t.constraints(p)
        if p ∈ ex_inds && !(Sc ⊆ S)
            @warn "Printing index $p as $S, but it is constrained to $Sc"
        end
    end

    translation = update_index_translation(t, translation)

    all_nonscalar_empty = isempty(t.sum_indices) && isempty(t.deltas) &&
                          isempty(t.tensors) && isempty(t.operators) &&
                          isempty(t.constraints)

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
            print_mo_index(io, t.constraints, translation, i)
        end
        print(io, '(')
        sep[] = false
    end

    for d in t.deltas
        printsep()
        print(io, (d, t.constraints, translation))
    end

    for ten in t.tensors
        printsep()
        print(io, (ten, t.constraints, translation))
    end

    for op in t.operators
        printsep()
        print(io, (op, t.constraints, translation))
    end

    if !isempty(t.sum_indices)
        print(io, ')')
    end

    if !do_index_translation
        printsep()

        print(io, "C(")

        isfirst = true

        for (i, s) in t.constraints
            if !isfirst
                print(io, ", ")
            end
            print_mo_index(io, t.constraints, translation, i)
            print(io, "∈", getshortname(s))
            isfirst = false
        end

        print(io, ')')
    end
end

# utility function to "copy" a term but replace the scalar with a new one
function new_scalar(t::Term{T1}, scalar::T2) where {T1<:Number,T2<:Number}
    Term(
        scalar,
        t.sum_indices,
        t.deltas,
        t.tensors,
        t.operators,
        t.constraints,
        t.max_simplified,
        true
    )
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
    operatortypes_a = [typeof(e) for e in a.operators]
    operatortypes_b = [typeof(e) for e in b.operators]

    tensorstrings_a = [get_symbol(t) for t in a.tensors]
    tensorstrings_b = [get_symbol(t) for t in b.tensors]

    (
        length(a.operators), operatortypes_a,
        length(a.sum_indices),
        length(a.tensors), tensorstrings_a,
        length(a.deltas),
        a.operators, a.sum_indices, a.tensors, a.deltas,
        a.constraints,
        -abs(a.scalar), -sign(a.scalar),
    ) < (
        length(b.operators), operatortypes_b,
        length(b.sum_indices),
        length(b.tensors), tensorstrings_b,
        length(b.deltas),
        b.operators, b.sum_indices, b.tensors, b.deltas,
        b.constraints,
        -abs(b.scalar), -sign(b.scalar),
    )
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

    t = copy(t, false)

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

function get_external_indices(t::Term)
    all_inds = get_all_indices(t)
    sort!(setdiff!(all_inds, t.sum_indices))
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

function possibly_equal_nonscalar(a::Term, b::Term)
    if (length(a.deltas), length(a.sum_indices)) !=
       (length(b.deltas), length(b.sum_indices))
        return false
    end

    if length(a.operators) != length(b.operators) ||
       !all(typeof(oa) == typeof(ob) for (oa, ob) in
            zip(a.operators, b.operators))
        return false
    end

    if length(a.tensors) != length(b.tensors) ||
       !all((get_symbol(ta), length(get_indices(ta))) ==
            (get_symbol(tb), length(get_indices(tb)))
            for (ta, tb) in zip(a.tensors, b.tensors))
        return false
    end

    true
end

function possibly_equal(a::Term, b::Term)
    if (a.scalar, length(a.deltas), length(a.sum_indices)) !=
       (b.scalar, length(b.deltas), length(b.sum_indices))
        return false
    end

    if length(a.operators) != length(b.operators) ||
       !all(typeof(oa) == typeof(ob) for (oa, ob) in
            zip(a.operators, b.operators))
        return false
    end

    if length(a.tensors) != length(b.tensors) ||
       !all((get_symbol(ta), length(get_indices(ta))) ==
            (get_symbol(tb), length(get_indices(tb)))
            for (ta, tb) in zip(a.tensors, b.tensors))
        return false
    end

    true
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

function constraints_could_be_equal_but_one(a::Term, b::Term)
    exinds = get_external_indices(a)

    if exinds != get_external_indices(b)
        return false
    end

    for p in exinds
        if a.constraints(p) != b.constraints(p)
            return false
        end
    end

    counts_a = Dict{Int,Int}()
    counts_b = Dict{Int,Int}()

    for (_, s) in a.constraints
        counts_a[s] = get(counts_a, s, 0) + 1
        counts_b[s] = 0
    end

    for (_, s) in b.constraints
        counts_b[s] = get(counts_b, s, 0) + 1
        if !haskey(counts_a, s)
            counts_a[s] = 0
        end
    end

    found_up = false
    found_down = false

    for (s, c) in counts_a
        Δ = c - get(counts_b, s, 0)
        if abs(Δ) >= 2
            return false
        elseif Δ == 1
            if found_up
                return false
            else
                found_up = true
            end
        elseif Δ == -1
            if found_down
                return false
            else
                found_down = true
            end
        end
    end

    found_up == found_down
end

function find_equal_perm(a::Term, b::Term)
    permuted_inds = copy(b.sum_indices)
    mapping = [p => p for p in permuted_inds]
    for perm in PermGen(length(b.sum_indices))
        copy!(permuted_inds, b.sum_indices)
        permute!(permuted_inds, perm.data)
        for (i, (p, q)) in enumerate(zip(b.sum_indices, permuted_inds))
            mapping[i] = p => q
        end

        new_b = sort_operators(exchange_indices(b, mapping))

        p = if equal_nonscalar(a, new_b)
            return true
        end

        p = if non_constraint_non_scalar_equal(a, new_b)
            constraints_equal_but_one(a, new_b)
        end

        if !isnothing(p)
            return (new_b, p)
        end
    end
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

    if s1 ⊊ s2
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
        min_t = min(min_t, sort_operators(exchange_indices(t, mapping)))
    end
    min_t
end

# This function is just a composition of other simplification functions
# in the recomended order to obtain a deterministic simplification of
# the term.
# TODO: This requires extensive testing of determinism and correctness.
# The specific steps and the order of them might need to be adjusted.
function simplify(t::Term)
    if t.max_simplified
        t
    else
        t |>
        lower_delta_indices |>
        simplify_summation_deltas |>
        lower_summation_indices |>
        sort_summation_indices |>
        sort_operators
    end
end

function set_max_simplified(t::Term)
    Term(
        t.scalar,
        t.sum_indices,
        t.deltas,
        t.tensors,
        t.operators,
        t.constraints,
        true,
        true
    )
end

function simplify_heavy(t::Term)
    if t.max_simplified
        t
    else
        t |>
        lower_delta_indices |>
        simplify_summation_deltas |>
        lower_summation_indices |>
        permute_all_sum_indices |>
        simplify_permute |>
        sort_operators |>
        set_max_simplified
    end
end

function split_indices(t::Term, (from, (to1, to2)))
    finished_terms = typeof(t)[]
    new_terms = typeof(t)[]
    old_terms = [t]

    while !isempty(old_terms)
        while !isempty(old_terms)
            t = pop!(old_terms)
            did_split = false
            for (p, s) in t.constraints
                if s == from
                    c1 = copy(t.constraints)
                    c2 = copy(t.constraints)

                    c1[p] = to1
                    c2[p] = to2

                    push!(new_terms, new_constraints(t, c1))
                    push!(new_terms, new_constraints(t, c2))
                    did_split = true
                    break
                end
            end
            if !did_split
                push!(finished_terms, t)
            end
        end

        old_terms, new_terms = new_terms, old_terms
    end

    finished_terms
end

# Some operator overloading (Not ment for external use):

function Base.:*(a::A, b::Term{B}) where {A<:Number,B<:Number}
    new_scalar(b, a * b.scalar)
end

# We are assuming commutative scalar multiplication
function Base.:*(a::Term{A}, b::B) where {A<:Number,B<:Number}
    b * a
end

# Returns 0 if constraints produce 0 and 1 otherwise
function fuse_constraints!(a::Constraints, b::Constraints)
    for (p, s) in b
        if haskey(a, p)
            if isdisjoint(a[p], s)
                return 0
            else
                a[p] = intersect(a[p], s)
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
        return Expression(zero(promote_type(A, B)))
    end

    b = make_space_for_indices(b, get_all_indices(a))
    a = make_space_for_indices(a, get_all_indices(b))

    Γ, e = reductive_commutator_fuse(a, b)

    if Γ == -1
        e - Expression([2 * b * a])
    else
        e
    end
end

function anticommutator(a::Term{A}, b::Term{B}) where {A<:Number,B<:Number}
    if isempty(a.operators) || isempty(b.operators)
        return Expression(zero(promote_type(A, B)))
    end

    b = make_space_for_indices(b, get_all_indices(a))
    a = make_space_for_indices(a, get_all_indices(b))

    Γ, e = reductive_commutator_fuse(a, b)

    if Γ == 1
        e + Expression([2 * b * a])
    else
        e
    end
end

function reductive_commutator(a::Term{A}, b::Term{B}) where
{A<:Number,B<:Number}
    if isempty(a.operators) || isempty(b.operators)
        return (1, Expression(zero(promote_type(A, B))))
    end

    b = make_space_for_indices(b, get_all_indices(a))
    a = make_space_for_indices(a, get_all_indices(b))

    reductive_commutator_fuse(a, b)
end

function reductive_commutator_fuse(a::Term{A}, b::Term{B}) where
{A<:Number,B<:Number}
    if isempty(a.operators) || isempty(b.operators)
        return (1, Expression(zero(promote_type(A, B))))
    end

    constraints = copy(a.constraints)
    if fuse_constraints!(constraints, b.constraints) == 0
        return (1, Expression(zero(promote_type(A, B))))
    end

    terms = Term{promote_type(A, B)}[]

    Γ = 1
    for j in eachindex(b.operators), i in reverse(eachindex(a.operators))
        δij, e = reductive_commutator(a.operators[i], b.operators[j])

        lhs = Operator[a.operators[1:i-1]; b.operators[1:j-1]]
        rhs = Operator[b.operators[j+1:end]; a.operators[i+1:end]]

        for t in e.terms
            constraints = copy(a.constraints)
            fuse_constraints!(constraints, t.constraints)
            fuse_constraints!(constraints, b.constraints)

            fused = Term(
                Γ * a.scalar * t.scalar * b.scalar,
                Int[a.sum_indices; t.sum_indices; b.sum_indices],
                KroneckerDelta[a.deltas; t.deltas; b.deltas],
                Tensor[a.tensors; t.tensors; b.tensors],
                Operator[lhs; t.operators; rhs],
                constraints
            )

            push!(terms, fused)

            Γ *= δij
        end
    end

    (Γ, Expression(terms))
end

# Function to express all operators in an expression in terms of
# elementary fermionic/bosinic anihilation and creation operators (if possible)
function convert_to_elementary_operators(t::Term{T}) where {T<:Number}
    ex = Expression([noop_part(t)])

    for o in t.operators
        ex = fuse(ex, convert_to_elementary_operators(o))
    end

    ex
end

function sort_operators(t::Term)
    t = copy(t)

    s = 1
    done = false
    while !done
        done = true
        for i in 1:length(t.operators)-1
            if t.operators[i] > t.operators[i+1]
                Γ, e = reductive_commutator(t.operators[i], t.operators[i+1])

                for i in eachindex(e.terms)
                    et = e[i]
                    e[i] = Term(
                        et.scalar,
                        et.sum_indices,
                        et.deltas,
                        et.tensors,
                        et.operators,
                        t.constraints
                    )
                end

                e = Expression(e.terms)

                if iszero(e)
                    done = false
                    s *= Γ
                    tmp = t.operators[i]
                    t.operators[i] = t.operators[i+1]
                    t.operators[i+1] = tmp
                end
            end
        end
    end

    s * t
end

export simplify_permute

# Function to do all permutations of the PermuteTensor and find the
# minimally-sorted term
# P_aibj F_bj -> P_aibj F_ai
function simplify_permute(t::Term)
    perm_tensors = [x for x in enumerate(t.tensors)
                    if x[2] isa PermuteTensor]
    if isempty(perm_tensors)
        return t
    end

    t = copy(t)

    sum_inds = copy(t.sum_indices)
    empty!(t.sum_indices)

    min_t = t
    for (p_i, tensor) in reverse!(perm_tensors)
        # Preallocate relevant arrays
        indices = get_indices(tensor)
        permuted_inds = copy(indices)
        mapping = [p => p for p in indices]
        N_pairs = length(indices) ÷ 2

        # Get all permutations of the pairs,
        # but drop the identity permutation.
        all_equal = true
        for perm in Iterators.drop(PermGen(N_pairs), 1)
            for (i, new_i) = enumerate(perm.data)
                permuted_inds[2new_i-1] = indices[2i-1]
                permuted_inds[2new_i] = indices[2i]
            end
            for (i, (p, q)) in enumerate(zip(indices, permuted_inds))
                mapping[i] = p => q
            end
            sorted_term = sort_operators(exchange_indices(t, mapping))
            if sorted_term != min_t
                all_equal = false
            end
            min_t = min(min_t, sorted_term)
        end

        if all_equal
            deleteat!(min_t.tensors, p_i)
            min_t = new_scalar(min_t, min_t.scalar * factorial(N_pairs))
        end
    end

    append!(min_t.sum_indices, sum_inds)

    return min_t
end

function Base.adjoint(t::Term)
    Term(
        t.scalar,
        t.sum_indices,
        t.deltas,
        t.tensors,
        [o' for o in Iterators.reverse(t.operators)],
        t.constraints,
        false
    )
end
