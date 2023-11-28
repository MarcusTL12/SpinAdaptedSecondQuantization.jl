export desymmetrize, make_permutation_mappings

function make_permutation_mappings(index_groups)
    mappings = Vector{Pair{Int,Int}}[]

    for perm in Iterators.drop(PermGen(length(index_groups)), 1)
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

function desymmetrize(ex::Expression{T}, mappings) where {T <: Number}
    accounted_for = Set{Int}()

    non_symmetric = Term{T}[]
    self_symmetric = Term{T}[]
    symmetrize = Term{T}[]

    for (i, t) in enumerate(ex.terms)
        if i ∈ accounted_for
            continue
        end

        other_inds = Int[]
        for mapping in mappings
            other_term = simplify_heavy(exchange_indices(t, mapping))
            if other_term == t
                push!(other_inds, i)
            else
                for (j, t2) in enumerate(ex.terms)
                    if j ∉ other_inds && other_term == t2
                        push!(other_inds, j)
                        break
                    end
                end
            end
        end

        if length(other_inds) < length(mappings)
            push!(non_symmetric, t)
            push!(accounted_for, i)
        elseif all(==(i), other_inds)
            push!(self_symmetric, t)
            push!(accounted_for, i)
        else
            push!(symmetrize, t)
            push!(accounted_for, i)
            union!(accounted_for, other_inds)
        end
    end

    non_symmetric = Expression(non_symmetric)
    self_symmetric = Expression(self_symmetric)
    symmetrize = Expression(symmetrize)

    symmetrize, self_symmetric, non_symmetric
end
