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

export scalar_type
function scalar_type(::Expression{T}) where {T<:Number}
    T
end

export summation, ∑, Σ

function summation(e::Expression, sum_indices)
    Expression([summation(t, sum_indices) for t in e.terms])
end

# Unicode alias
∑(e, s) = summation(e, s)
Σ(e, s) = summation(e, s)

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

"""
    simplify(ex::Expression)

Applies various simplification strategies to the given expression.
The most important are:
- Remove Kronecker deltas when indices are summed over.
- Perform various operations to try to lower the lexiographic ordering of the
    each term including:
    - Lowering indices showing up in a Kronecker delta to the lowest of the
        indices in said delta.
    - Lowering summation indices to lowest available indices.
    - Reorder summation indices to lower the lexiographic ordering of the order
        they show up within the term.
    - Bubble sort operators swaping neighbouring operators if:
        - They commute
        - The swap would lower the terms lexiographic ordering.
- Look for pairs of terms differing only in index constraints and scalar
    prefactors, and tries to combine them.

# Examples:

Collapsing deltas:
```jldoctest simplify
julia> using SpinAdaptedSecondQuantization

julia> ∑(real_tensor("h", 1, 2) * electron(1, 2) * δ(1, 2), [2])
∑_q(δ_pq h_pq)
julia> simplify(ans)
h_pp
```

Combining terms differing only by index constraints and scalars:
```jldoctest simplify
julia> real_tensor("h", 1, 2) * (occupied(1, 2) + occupied(1) * virtual(2))
h_ia
+ h_ij
julia> simplify(ans)
h_ip
```
```jldoctest simplify
julia> ∑(real_tensor("h", 1, 2) * E(1, 2) * electron(1, 2), [1, 2])
∑_pq(h_pq E_pq)
julia> ans - ∑(real_tensor("h", 1, 2) * E(1, 2) * electron(1) * occupied(2), [1, 2])
∑_pq(h_pq E_pq)
- ∑_pi(h_pi E_pi)
julia> simplify(ans)
∑_pa(h_pa E_pa)
```

Sorting commuting operators:
```jldoctest simplify
julia> E(3, 4) * E(1, 2) * electron(1:4...)
E_rs E_pq
julia> simplify(ans) # not able to swap as commutator is not zero
E_rs E_pq
julia> ans * occupied(2, 4) * virtual(1, 3)
E_bj E_ai
julia> simplify(ans) # now able to swap
E_ai E_bj
```
```jldoctest simplify
julia> E(3, 4) * E(5, 6) * E(1, 2) * electron(5, 6) * occupied(2, 4) * virtual(1, 3)
E_bj E_pq E_ai
julia> simplify(ans) # Can not swap first and last because of non commuting between
E_bj E_pq E_ai
julia> E(5, 6) * E(3, 4) * E(1, 2) * electron(5, 6) * occupied(2, 4) * virtual(1, 3)
E_pq E_bj E_ai
julia> simplify(ans) # Can swap adjacent last two, but not first
E_pq E_ai E_bj
```
"""
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

export simplify_heavy

"""
    simplify_heavy(ex::Expression,
        [mapping=GeneralOrbital => (OccupiedOrbital, VirtualOrbital)])

Similar to [`simplify`](@ref), but enabling a few extra simplification
strategies that can be quite expensive, especially as the number of summation
indices in a term grows large.

The most expensive step is to iterate over *every permutation* of
the summation indices of each term in order to find the order that has the least
lexiographic ordering, which naturally has a scaling of O(n!) where n is the
larget number of summation indices in a single term. This is quite managable for
"small" numbers of indices (less than ~8), but quickly becomes unfeasible.

!!! note
    The iteration over permutations of summation indices happens after summation
    indices that show up in Kronecker deltas, so if your terms contain many
    Kronecker deltas this might yet be possible to run.

The optional `mapping` argument defaults to splitting terms with general indices
into a term with the index being occupied plus a term where the index is
virtual, then performing simplification when all terms are split and recombining
terms according to the same splitting if possible. This often allows for a few
additional cancellations to occur than using only [`simplify`](@ref)
If dealing with other index spaces, this splitting can be specified can
be specified to be some other set of index spaces. See [`split_indices`](@ref).

# Examples:

HF energy with full real orbital symmetry:
```jldoctest simplify_heavy
julia> using SpinAdaptedSecondQuantization

julia> h = ∑(rsym_tensor("h", 1, 2) * E(1, 2) * electron(1, 2), 1:2)
∑_pq(h_pq E_pq)
julia> g = 1//2 * \
∑(rsym_tensor("g", 1:4...) * e(1:4...) * electron(1:4...), 1:4)
1/2 ∑_pqrs(g_pqrs E_pq E_rs)
- 1/2 ∑_pqrs(δ_qr g_pqrs E_ps)
julia> H = simplify(h + g + real_tensor("h_nuc"));

julia> E_HF = simplify(hf_expectation_value(H)) # Not nice
h_nuc
+ 2 ∑_i(h_ii)
+ 2 ∑_ij(g_iijj)
- ∑_pi(g_pipi)
+ ∑_ia(g_iaia)
julia> simplify_heavy(E_HF) # Much nicer!
h_nuc
+ 2 ∑_i(h_ii)
+ 2 ∑_ij(g_iijj)
- ∑_ij(g_ijij)
```

Simplification by combining terms differing by index constraints
```jldoctest simplify_heavy
julia> ∑(real_tensor("h", 1, 2) * E(1, 2) * electron(1, 2), [1, 2])
∑_pq(h_pq E_pq)
julia> ans - ∑(real_tensor("h", 1, 2) * E(1, 2) * occupied(1, 2), [1, 2])
∑_pq(h_pq E_pq)
- ∑_ij(h_ij E_ij)
julia> ans - ∑(real_tensor("h", 1, 2) * E(1, 2) * virtual(1, 2), [1, 2])
∑_pq(h_pq E_pq)
- ∑_ab(h_ab E_ab)
- ∑_ij(h_ij E_ij)
julia> simplify(ans)
∑_pq(h_pq E_pq)
- ∑_ab(h_ab E_ab)
- ∑_ij(h_ij E_ij)
julia> simplify_heavy(ans)
∑_ai(h_ai E_ai)
+ ∑_ia(h_ia E_ia)
```
"""
function simplify_heavy(ex::Expression,
    mapping=GeneralOrbital => (OccupiedOrbital, VirtualOrbital))
    done = false
    ex = simplify(ex)
    while !done
        new_ex = split_indices(ex, mapping) |>
                 simplify_heavy_terms |>
                 try_add_constraints
        done = new_ex == ex
        ex = new_ex
    end
    ex
end

export sort_operators
function sort_operators(ex::Expression)
    Expression([sort_operators(t) for t in ex.terms])
end

"""
    split_indices(ex::Expression, mapping)

Splits the indices of each term according to the given mapping.

```jldoctest @example
julia> using SpinAdaptedSecondQuantization

julia> h = ∑(real_tensor("h", 1, 2) * E(1, 2) * electron(1, 2), 1:2)
∑_pq(h_pq E_pq)

julia> SASQ.split_indices(h,\
 GeneralOrbital => (OccupiedOrbital, VirtualOrbital))
∑_ab(h_ab E_ab)
+ ∑_ai(h_ai E_ai)
+ ∑_ia(h_ia E_ia)
+ ∑_ij(h_ij E_ij)
```
"""
function split_indices(ex::Expression, mapping)
    terms = eltype(ex.terms)[]

    for t in ex.terms
        append!(terms, split_indices(t, mapping))
    end

    Expression(terms)
end

export commutator

"""
    commutator(a::Expression, b::Expression)

Compute commutator between `a` and `b`. This is computed termwise using the
[`reductive_commutator`](@ref) function, then adding/subtracting the swapped
result if the reductive commutator for a given pair of terms turned out to be
an anti-commutator.
"""
function commutator(a::Expression{A}, b::Expression{B}) where
{A<:Number,B<:Number}
    terms = Term{promote_type(A, B)}[]

    for t1 in a.terms, t2 in b.terms
        append!(terms, commutator(t1, t2).terms)
    end

    Expression(terms)
end

export anticommutator

"""
    anticommutator(a::Expression, b::Expression)

Compute commutator between `a` and `b`. This is computed termwise using the
[`reductive_commutator`](@ref) function, then adding/subtracting the swapped
result if the reductive commutator for a given pair of terms turned out to be
an commutator.
"""
function anticommutator(a::Expression{A}, b::Expression{B}) where
{A<:Number,B<:Number}
    terms = Term{promote_type(A, B)}[]

    for t1 in a.terms, t2 in b.terms
        append!(terms, anticommutator(t1, t2).terms)
    end

    Expression(terms)
end

"""
    commutator(A::Expression, B::Expression, n::Integer)

Compute `n` nested commutators ``[...[[[A], B], B]..., B]``.
"""
function commutator(A::Expression, B::Expression, n::Integer)
    # [A, B]_n
    # [A, B]_3 = [[[A, B], B], B]
    X = A
    for _ = 1:n
        X = commutator(X, B)
    end
    return X
end

"""
    commutator(A, B, C, ...)

Computes the nested commutator ``[...[[A, B], C], ...]``
"""
function commutator(A, B, rest...)
    commutator(commutator(A, B), rest...)
end

"""
    commutator(A, Bs)

Compute nested commutator with all expressions in iteratable Bs.

``
[...[[[A], B_1], B_2], ...]
``
"""
function nested_commutator(A, Bs)
    for B in Bs
        A = commutator(A, B)
    end
    A
end

function bch_smart_kernel(A, Bs, n)
    if n == 0
        return A
    end

    inds = zeros(Int, n)

    acc = Term[]

    function calculate_prefactor()
        cur_group_start = 1
        cur_group_ind = inds[1]
        fac = 1
        for i in 2:n
            if inds[i] != cur_group_ind
                fac *= factorial(i - cur_group_start)
                cur_group_start = i
                cur_group_ind = inds[i]
            end
        end

        fac *= factorial(n + 1 - cur_group_start)

        fac
    end

    function rec(i)
        if i == n + 1
            choices = [Bs[j] for j in inds]
            c = nested_commutator(A, choices)
            append!(acc, (1 // calculate_prefactor() * c).terms)
            return
        end

        last_max = if i == 1
            length(Bs)
        else
            inds[i-1]
        end

        for j in 1:last_max
            inds[i] = j
            rec(i + 1)
        end
    end

    rec(1)

    Expression(acc)
end

"""
    bch(A, Bs::AbstractArray, n)

Equivalent to `bch(A, sum(Bs), n)` if all `Bs` commute with eachother, but
evalutating in a smart way in order to avoid computing equivalent terms many
times, for example instead of computing both ``\frac12 [[A, B_1], B_2]`` and
``\frac12 [[A, B_2], B_1]``, it will only compute the first one and multiply
by 2. This can be much faster if one has many commuting terms in `Bs`.

!!! warning
    Should not be used if some expressions in `Bs` do not commute.
"""
function bch(A, Bs::AbstractArray, n)
    acc = Term[]

    for i in 0:n
        append!(acc, bch_smart_kernel(A, Bs, i).terms)
    end

    Expression(acc)
end

"""
    bch(A, B, n)

Compute Baker-Campbell-Haussdorff expansion

``
e^{-B} A e^B
=
A + [A, B] + \\frac{1}{2!} [[A, B], B] + \\frac{1}{3!} [[[A, B], B], B] +
\\frac{1}{n!} [...[[[A, B], B], B], ...B]
``

up to nth order. Terminates if commutator becomes zero.
"""
function bch(A, B::Expression, n)
    # Baker-Campbell-Haussdorff expansion,
    # e^-B A e^B = A + 1/1! [A,B] + 1/2! [[A,B],B] + ... + 1/n! [A,B]_n
    X = A
    Y = A
    for i = 1:n
        Y = 1 // i * commutator(Y, B)
        if iszero(Y)
            break
        end
        X += Y
    end
    X
end

export act_on_bra

"""
    act_on_bra(x::Expression, [max_ops=Inf])

Projects `x` on a HF bra by projecting the adjoint on a ket and ajointing the
result.

```julia
act_on_bra(x, max_ops=Inf) = act_on_ket(x', max_ops)'
```

See [`act_on_ket`](@ref)
"""
act_on_bra(x, max_ops=Inf) = act_on_ket(x', max_ops)'

export act_eT_on_bra

"""
    act_eT_on_bra(braop::Expression, T::Expression; [max_n=Inf], [max_ops=Inf])

computes

``
\\langle\\text{HF}| \\text{braop}\\ e^T
``

This can be useful for evaluating certain coupled cluster expressions from
left to right. Since the cluster operator T is purely exciting acting ``e^T``
to the left will terminate after only a few terms depending on the left state
being projected on.

!!! warning
    If `T` contains de-excitations, the expression will not truncate and the
    parameter `max_n` must be specified with a finite positive integer for
    the function to terminate, truncating the expression.
    `max_n` specifies the highest order term from the
    taylor expansion of `e^T` to include. If `T` is purely exciting, the
    expansion exits early when it terminates.
"""
function act_eT_on_bra(braop, T; max_n=Inf, max_ops=Inf)
    i = 0
    acc = act_on_bra(braop)
    proj = acc

    while i <= max_n
        i += 1
        proj = simplify(act_on_bra(proj * T)) * 1 // i
        if iszero(proj)
            break
        end
        acc += proj
    end

    simplify(act_on_bra(acc, max_ops))
end

export convert_to_elementary_operators

"""
    convert_to_elementary_operators(ex::Expression)

converts all operators to their representation using "elementary"
annihilation and creation operators if implemented for all operator types in
`ex`.
"""
function convert_to_elementary_operators(ex::Expression{T}) where {T<:Number}
    terms = Term{T}[]

    for t in ex.terms
        append!(terms, convert_to_elementary_operators(t).terms)
    end

    Expression(terms)
end

"""
    Base.adjoint(ex::Expression)

Termwise adjoint ex. Equivalent to the postfix notation `ex'`

# Examples

```jldoctest
julia> using SpinAdaptedSecondQuantization

julia> Eai(a, i) = E(a, i) * occupied(i) * virtual(a)
Eai (generic function with 1 method)
julia> T2 = 1//2 * ∑(psym_tensor("t", 1:4...) * Eai(1, 2) * Eai(3, 4), 1:4)
1/2 ∑_aibj(t_aibj E_ai E_bj)
julia> adjoint(T2)
1/2 ∑_aibj(t_aibj E_jb E_ia)
julia> T2'
1/2 ∑_aibj(t_aibj E_jb E_ia)
julia> simplify(ans)
1/2 ∑_iajb(t_aibj E_ia E_jb)
```
"""
function Base.adjoint(ex::Expression)
    Expression([t' for t in ex.terms])
end

"""
    Base.:(^)(ex::Expression, n::Integer)

Computes `ex^n` by repeated multiplication.

!!! warning
    n must be a non negative integer

# Examples

```jldoctest
julia> using SpinAdaptedSecondQuantization

julia> h = ∑(real_tensor("h", 1, 2) * E(1, 2) * electron(1, 2), 1:2)
∑_pq(h_pq E_pq)
julia> h^2
∑_pqrs(h_pq h_rs E_pq E_rs)
julia> h^3
∑_pqrstu(h_pq h_rs h_tu E_pq E_rs E_tu)
```
"""
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
