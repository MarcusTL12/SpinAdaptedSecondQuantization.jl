export look_for_tensor_replacements, make_exchange_transformer

function do_tensor_replacement(t::Term, transformer)
    old_tensors = t.tensors

    new_terms = Term[]
    other_terms = Term[]

    for (i, tens) in enumerate(old_tensors)
        result = transformer(tens)

        if !isnothing(result)
            new_scal, new_tens, other_scal, other_tens = result

            new_tensors = copy(old_tensors)
            new_tensors[i] = new_tens

            new_term = Term(
                t.scalar * new_scal,
                t.sum_indices,
                t.deltas,
                new_tensors,
                t.operators,
                t.constraints
            )

            other_tensors = copy(old_tensors)
            other_tensors[i] = other_tens

            other_term = Term(
                t.scalar * other_scal,
                t.sum_indices,
                t.deltas,
                other_tensors,
                t.operators,
                t.constraints
            )

            push!(new_terms, new_term)
            push!(other_terms, other_term)
        end
    end

    new_terms, other_terms
end

function make_exchange_transformer(from, to)
    function g2L_transformer(t::T) where {T<:Tensor}
        if length(get_indices(t)) != 4 || get_symbol(t) != from
            return
        end

        -1 // 2, reorder_indices(t, [1, 4, 3, 2]), 1 // 2, T(to, get_indices(t))
    end

    g2L_transformer
end

function look_for_tensor_replacements(ex::Expression, transformer)
    is_done = false

    while !is_done
        is_done = true

        for i in eachindex(ex.terms)
            replacements, other_replacements =
                do_tensor_replacement(ex[i], transformer)

            for (replacement, other_replacement) in
                zip(replacements, other_replacements)
                for j in eachindex(ex.terms)
                    if i != j
                        if possibly_equal(ex[j], replacement)
                            simple_replacement = simplify_heavy(replacement)
                            if ex[j] == simple_replacement
                                is_done = false

                                new_terms = copy(ex.terms)
                                new_terms[i] = simplify_heavy(other_replacement)
                                deleteat!(new_terms, j)
                                ex = Expression(new_terms)
                            end
                        end
                    end
                    if !is_done
                        break
                    end
                end
                if !is_done
                    break
                end
            end
            if !is_done
                break
            end
        end
    end

    ex
end
