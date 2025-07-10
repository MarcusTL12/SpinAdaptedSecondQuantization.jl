export pick_biorthogonal

function project_biorthogonal_operator(
    o::SingletExcitationOperator,
    template::SingletExcitationOperator)

    [
        KroneckerDelta(o.p, template.p),
        KroneckerDelta(o.q, template.q),
    ]
end

function pick_biorthogonal(
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
