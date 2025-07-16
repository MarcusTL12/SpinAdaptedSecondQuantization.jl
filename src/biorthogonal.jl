export project_biorthogonal

function project_biorthogonal_operator(
    o::SingletExcitationOperator,
    template::SingletExcitationOperator)

    [
        KroneckerDelta(o.p, template.p),
        KroneckerDelta(o.q, template.q),
    ]
end

"""
    project_biorthogonal(ex::Expression, template::Expression)

This has the effect of projecting the expression `ex` on the biorthogonal bra
state of the `template` expression.

For a singlet double excitation, the following will be equivalent:

```
bra = (
        1//3 * E(1, 2) * E(3, 4) +
        1//6 * E(1, 4) * E(3, 2)
    )' * occupied(2, 4) * virtual(1, 3)

result = hf_expectation_value(bra * ex)
```

```
ket = E(1, 2) * E(3, 4) * occupied(2, 4) * virtual(1, 3)
result = project_biorthogonal(act_on_ket(ex), ket)
result = symmetrize(result, make_permutation_mappings([(1, 2), (3, 4)]))
```
"""
function project_biorthogonal(
    x::Expression{T}, template_ex::Expression) where {T<:Number}
    @assert length(template_ex.terms) == 1

    template = template_ex[1]

    @assert isempty(template.deltas)
    @assert isempty(template.tensors)
    @assert isempty(template.sum_indices)
    @assert isone(template.scalar)

    out_inds = get_all_indices(template)

    template_types = [typeof(o) for o in template.operators]

    terms = Term{T}[]

    for t in x.terms
        if length(t.operators) == length(template_types) &&
           all(typeof(o) == tp for (o, tp) in zip(t.operators, template_types))
            t = make_space_for_indices(copy(t), out_inds)
            for (o, o_template) in zip(t.operators, template.operators)
                new_deltas = project_biorthogonal_operator(o, o_template)

                append!(t.deltas, new_deltas)
            end

            empty!(t.operators)

            push!(terms, t)
        end
    end

    simplify_heavy(Expression(terms))
end
