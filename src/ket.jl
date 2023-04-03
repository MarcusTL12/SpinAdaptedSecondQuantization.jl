export act_on_ket

function act_on_ket(ex::Expression{T}, max_ops=Inf) where {T}
    nth = Threads.nthreads()
    terms = [Term{T}[] for _ in 1:nth]
    Threads.@threads for id in 1:nth
        for i in id:nth:length(ex.terms)
            append!(terms[id], act_on_ket(ex[i], max_ops).terms)
        end
    end

    all_terms, rest = Iterators.peel(terms)
    for other_terms in rest
        append!(all_terms, other_terms)
    end

    Expression(all_terms)
end

function act_on_ket(t::Term{A}, max_ops) where {A<:Number}
    if iszero(t.scalar)
        return Expression(zero(A))
    end
    if isempty(t.operators)
        return Expression([t])
    end

    copyt = copy(t)
    right_op = pop!(copyt.operators)
    right_op_act = act_on_ket(right_op)
    copyt_act = act_on_ket(copyt,
        max_ops - minimum(length(t.operators) for t in right_op_act.terms))

    terms = Term{A}[]
    for r in right_op_act.terms
        Γ, comm = reductive_commutator_fuse(copyt, r)

        if length(r.operators) <= max_ops
            new_max = max_ops - length(r.operators)
            append!(terms, Γ * fuse(r, ter)
                           for ter in copyt_act.terms
                           if length(ter.operators) <= new_max)
        end

        append!(terms, act_on_ket(comm, max_ops).terms)
    end

    Expression(terms)
end

export hf_expectation_value
hf_expectation_value(ex::Expression) = act_on_ket(ex, 0)
