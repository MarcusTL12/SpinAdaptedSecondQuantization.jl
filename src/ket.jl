export act_on_ket

function act_on_ket(ex :: Expression)
    sum(act_on_ket(t) for t in ex.terms)
end

function act_on_ket(t :: Term)
    if iszero(t.scalar)
        return Expression(0)
    end
    if isempty(t.operators)
        return Expression([t])
    end

    copyt = copy(t)
    right_op = pop!(copyt.operators)
    right_op_act = act_on_ket(right_op)
    copyt_act = act_on_ket(copyt)

    ex = Expression(0)
    for r in right_op_act.terms
        ex += Expression([fuse(r, ter) for ter in copyt_act.terms])
        ex += act_on_ket(commutator_fuse(copyt, r))
    end

    return ex
end

function act_on_ket(op :: SingletExcitationOperator)
    p = op.p
    q = op.q
    E(p, q) * virtual(p) * occupied(q) +
        2 * δ(p, q) * occupied(p, q)
end