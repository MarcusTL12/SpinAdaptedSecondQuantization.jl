export hf_expectation_value

function hf_expectation_value(ex::Expression{T}) where {T<:Number}
    terms = Term{T}[]

    for t in ex.terms
        append!(terms, hf_expectation_value(t).terms)
    end

    Expression(terms)
end

function hf_expectation_value(t::Term)
    op_ex = hf_expectation_value(t.operators, t.constraints)

    Expression([fuse(noop_part(t), t2) for t2 in op_ex.terms])
end

function hf_expectation_value(ops::AbstractVector{Operator},
    constraints::Constraints)
    if isempty(ops)
        return Expression(1)
    elseif length(ops) == 1
        return hf_expectation_value(only(ops))
    end

    # Dynamic dispatch on type of first operator
    hf_expectation_value(first(ops), (@view ops[2:end]), constraints)
end

# Expectation value trick for specific first operator
function hf_expectation_value(firstop::SingletExcitationOperator, rest,
    constraints)
    if constraints(firstop.p) <: VirtualOrbital
        return zero(Expression{Int64})
    else
        if constraints(firstop.q) <: OccupiedOrbital
            hf_expectation_value(firstop) *
            hf_expectation_value(rest, constraints)
        elseif constraints(firstop.q) <: VirtualOrbital
            hf_expectation_value(
                commutator(
                    Expression(firstop),
                    Expression(rest)
                ) * constrain(constraints)
            )
        else
            assume_occ = hf_expectation_value(firstop) *
                         hf_expectation_value(rest, constraints) *
                         occupied(firstop.q)

            assume_vir = hf_expectation_value(commutator(
                Expression(firstop) *
                virtual(firstop.q),
                Expression(Vector(rest))
            ))

            assume_occ + assume_vir
        end
    end
end

function hf_expectation_value(o::SingletExcitationOperator)
    2δ(o.p, o.q) * occupied(o.p, o.q)
end
