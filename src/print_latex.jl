export latex

function latex(ex::Expression)
    io = IOBuffer()

    if isempty(ex.terms)
        throw("Expressions should not be empty,\
but rather include a single zero term")
    end

    t, rest = Iterators.peel(ex.terms)

    if t.scalar < 0
        print_latex(io, "- ", new_scalar(t, -t.scalar))
    else
        print_latex(io, t)
    end

    for t in rest
        if t.scalar < 0
            print_latex(io, "\n- ", new_scalar(t, -t.scalar))
        else
            print_latex(io, "\n+ ", t)
        end
    end

    String(take!(io))
end

function print_latex(io::IO, t::Term)

end
