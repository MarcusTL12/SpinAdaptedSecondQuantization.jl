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

    newt = copy(t)
    right_op = pop!(newt.operators)
    right_op_acted = act_on_ket(right_op)
    newt_act = act_on_ket(newt)

    t1 = Expression([fuse(right_op_acted[1], ter) for ter in newt_act.terms])
    t2 = Expression([fuse(right_op_acted[2], ter) for ter in newt_act.terms])
    t3 = act_on_ket(commutator(newt, right_op_acted[1]))
    t4 = act_on_ket(commutator(newt, right_op_acted[2]))

    return t1 + t2 + t3 + t4
end

function act_on_ket(op :: SingletExcitationOperator)
    p = op.p
    q = op.q
    [(E(p, q) * virtual(p) * occupied(q))[1],
     (2 * δ(p, q) * occupied(p, q))[1]]
end