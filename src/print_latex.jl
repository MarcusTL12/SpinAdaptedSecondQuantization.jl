export latex

latex(ex::Expression) = latex(ex, IndexTranslation())

function latex(ex::Expression, translation::IndexTranslation)
    io = IOBuffer()

    if isempty(ex.terms)
        throw("Expressions should not be empty,\
but rather include a single zero term")
    end

    t, rest = Iterators.peel(ex.terms)

    if t.scalar < 0
        print(io, "- ")
        print_latex(io, new_scalar(t, -t.scalar), translation)
    else
        print_latex(io, t, translation)
    end

    for t in rest
        if t.scalar < 0
            print(io, "\n- ")
            print_latex(io, new_scalar(t, -t.scalar), translation)
        else
            print(io, "\n+ ")
            print_latex(io, t, translation)
        end
    end

    String(take!(io))
end

function print_latex(io::IO, t::Term, translation::IndexTranslation)
    sep = Ref(false)

    function printsep()
        if sep[]
            print(io, ' ')
        end
        sep[] = true
    end

    ex_inds = get_external_indices(t)

    for (p, (S, _)) in translation
        Sc = t.constraints(p)
        if p ∈ ex_inds && !(Sc <: S)
            @warn "Printing index $p as $S, but it is constrained to $Sc"
        end
    end

    translation = update_index_translation(t, translation)

    all_nonscalar_empty = isempty(t.sum_indices) && isempty(t.deltas) &&
                          isempty(t.tensors) && isempty(t.operators) &&
                          isempty(t.constraints)

    if !isone(t.scalar)
        if isone(-t.scalar)
            print(io, '-')
        else
            printsep()
            printscalar(io, t.scalar)
        end
    elseif all_nonscalar_empty
        printscalar(io, t.scalar)
        sep[] = true
    end

    if !isempty(t.sum_indices)
        printsep()
        print(io, "\\sum_{")
        for i in t.sum_indices
            print_latex_mo_index(io, t.constraints, translation, i)
        end
        print(io, "}{")
        sep[] = false
    end

    for d in t.deltas
        printsep()
        print_latex(io, (d, t.constraints, translation))
    end

    for ten in t.tensors
        printsep()
        print_latex(io, (ten, t.constraints, translation))
    end

    for op in t.operators
        printsep()
        print_latex(io, (op, t.constraints, translation))
    end

    if !isempty(t.sum_indices)
        print(io, '}')
    end

    constraint_noprint = index_color_latex ? get_non_constraint_indices(t) : Int[]
    filter!(constraint_noprint) do x
        haskey(colors, t.constraints(x))
    end
    constraint_print = [i for (i, _) in t.constraints if i ∉ constraint_noprint]
    filter!(constraint_print) do x
        is_strict_subspace(t.constraints(x), translation(x)[1])
    end

    if !isempty(constraint_print)
        printsep()

        print(io, "C(")

        isfirst = true

        for i in constraint_print
            if !isfirst
                print(io, ", ")
            end
            print_latex_mo_index(io, t.constraints, translation, i)
            print(io, "\\in ", getshortname(t.constraints(i)))
            isfirst = false
        end

        print(io, ')')
    end
end
