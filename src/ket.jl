export act_on_ket

function act_on_ket(ex::Expression{T}) where {T}
    nth = Threads.nthreads()
    terms = [Term{T}[] for _ in 1:nth]
    Threads.@threads for id in 1:nth
        for i in id:nth:length(ex.terms)
            append!(terms[id], act_on_ket(ex[i]).terms)
        end
    end

    all_terms, rest = Iterators.peel(terms)
    for other_terms in rest
        append!(all_terms, other_terms)
    end

    Expression(all_terms)
end

function act_on_ket(t::Term{A}) where {A<:Number}
    if iszero(t.scalar)
        return Expression(zero(A))
    end
    if isempty(t.operators)
        return Expression([t])
    end

    copyt = copy(t)
    right_op = pop!(copyt.operators)
    right_op_act = act_on_ket(right_op)
    copyt_act = act_on_ket(copyt)

    terms = Term{A}[]
    for r in right_op_act.terms
        append!(terms, fuse(r, ter) for ter in copyt_act.terms)
        append!(terms, act_on_ket(commutator_fuse(copyt, r)).terms)
    end

    Expression(terms)
end

function act_on_ket(op::SingletExcitationOperator)
    p = op.p
    q = op.q
    E(p, q) * virtual(p) * occupied(q) +
    2 * δ(p, q) * occupied(p, q)
end