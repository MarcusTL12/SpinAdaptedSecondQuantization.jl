export symmetrize, desymmetrize, make_permutation_mappings

function make_permutation_mappings(index_groups)
    mappings = Vector{Pair{Int,Int}}[]

    for perm in PermGen(length(index_groups))
        pd = perm.data
        mapping = Pair{Int,Int}[]

        for (i, j) in enumerate(pd)
            for (from, to) in zip(index_groups[i], index_groups[j])
                push!(mapping, from => to)
            end
        end

        push!(mappings, mapping)
    end

    mappings
end

function desymmetrize(ex::Expression{T}, mappings) where {T<:Number}
    accounted_for = Set{Int}()

    non_symmetric = Term{T}[]
    self_symmetric = Term{T}[]
    symmetrize = Term{T}[]

    for (i, t) in enumerate(ex.terms)
        if i ∈ accounted_for
            continue
        end

        other_inds = Dict{Int,Int}()
        for mapping in mappings
            other_term = simplify_heavy(exchange_indices(t, mapping))
            if other_term == t
                other_inds[i] = get(other_inds, i, 0) + 1
            else
                for (j, t2) in enumerate(ex.terms)
                    if other_term == t2
                        other_inds[j] = get(other_inds, j, 0) + 1
                        break
                    end
                end
            end
        end

        if all(==(i), keys(other_inds))
            push!(self_symmetric, t)
            push!(accounted_for, i)
        elseif allequal(values(other_inds))
            push!(symmetrize, t * (1 // first(values(other_inds))))
            push!(accounted_for, i)
            union!(accounted_for, keys(other_inds))
        else
            push!(non_symmetric, t)
            push!(accounted_for, i)
        end
    end

    non_symmetric = Expression(non_symmetric)
    self_symmetric = Expression(self_symmetric)
    symmetrize = Expression(symmetrize)

    symmetrize, self_symmetric, non_symmetric
end

function symmetrize(ex::Expression{T}, mappings) where {T<:Number}
    terms = Term{T}[]

    for t in ex.terms
        for mapping in mappings
            other_term = exchange_indices(t, mapping)
            push!(terms, other_term)
        end
    end

    Expression(terms)
end
