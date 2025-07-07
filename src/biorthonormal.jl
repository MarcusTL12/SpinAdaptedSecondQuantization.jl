export pick_biorthonormal

function project_biorthonormal_term(t)
    t = copy(t)

    n_op = length(t.operators)

    out_inds = 1:2n_op

    t = SASQ.make_space_for_indices(t, out_inds)

    for (base_ind, o) in enumerate(t.operators)
        if !(o isa SingletExcitationOperator)
            throw("Biorthonormal currently only " *
                  "supported for Epq type operators only")
        end

        a = 2base_ind - 1
        i = 2base_ind

        push!(t.deltas, SASQ.KroneckerDelta(a, o.p))
        push!(t.deltas, SASQ.KroneckerDelta(i, o.q))
    end

    empty!(t.operators)

    t
end

function pick_biorthonormal(x)
    tiers = Dict()

    for t in x.terms
        n = length(t.operators)

        tiers[n] = get(tiers, n, SASQ.Expression(0)) +
                   SASQ.Expression([project_biorthonormal_term(t)])
    end

    [simplify_heavy(tiers[n-1]) for n in 2:length(tiers)]
end
