export print_code, print_eT_code, print_latex
using Printf

""" Print term to numpy's einsum routine
"""
function print_code(t::Term, symbol::String, translation)
    # Print python np.einsum code
    scalar_str = @sprintf "%+12.8f" t.scalar
    translation = update_index_translation(t, translation)

    external_int = get_external_indices(t)
    external = sprint(SASQ.print_mo_index, t.constraints, translation, external_int...)

    temp_color = index_color
    disable_color()

    # Make einsum_str
    einsum_str = "\""
    for a in t.tensors
        einsum_str *= sprint(SASQ.print_mo_index, t.constraints, translation, get_indices(a)...)
        einsum_str *= ","
    end
    einsum_str = einsum_str[begin:end-1] * "->" * external * "\""

    if (temp_color)
        enable_color()
    end

    # Make tensor_str
    tensor_str = ""
    for a in t.tensors
        tensor_str *= ", $(get_symbol(a))_"
        for b in get_indices(a)
            if t.constraints[b] == VirtualOrbital
                tensor_str *= "v"
            elseif t.constraints[b] == OccupiedOrbital
                tensor_str *= "o"
            else
                throw("Space not supported")
            end
        end
    end

    "$(symbol)_$(external) += $scalar_str * np.einsum($einsum_str$tensor_str, optimize=\"optimal\");"
end

""" Used with https://github.com/alexancp/einsumpath-to-eT to generate Fortran subroutines for eT.
"""
function print_eT_code(t::Term, symbol, translation, routine_name)
    prefactor = @sprintf "%+12.8f" t.scalar
    translation = update_index_translation(t, translation)

    external_int = get_external_indices(t)
    external = sprint(SASQ.print_mo_index, t.constraints, translation, external_int...)

    temp_color = index_color
    disable_color()

    # Make einsum_str
    einsum_str = ""
    for a in t.tensors
        einsum_str *= sprint(SASQ.print_mo_index, t.constraints, translation, get_indices(a)...)
        einsum_str *= ","
    end
    einsum_str = einsum_str[begin:end-1] * "->" * external

    # Make tensor_str
    tensor_str = ""
    for a in t.tensors
        tensor_str *= ", $(get_symbol(a))_"
        for b in get_indices(a)
            if t.constraints[b] == VirtualOrbital
                tensor_str *= "v"
            elseif t.constraints[b] == OccupiedOrbital
                tensor_str *= "o"
            else
                throw("Space not supported")
            end
        end
    end

    # Make symbol_str
    symbol_str = ""
    for a in t.tensors
        symbol_str *= ", \"$(get_symbol(a))_"
        for b in get_indices(a)
            if t.constraints[b] == VirtualOrbital
                symbol_str *= "v"
            elseif t.constraints[b] == OccupiedOrbital
                symbol_str *= "o"
            else
                throw("Space not supported")
            end
        end
        symbol_str *= '\"'
    end
    tensor_str = tensor_str[3:end]
    symbol_str = symbol_str[3:end]

    symbol_str = replace(symbol_str, "t_vovovo" => "t3", "t_vovo" => "t2")

    if (temp_color)
        enable_color()
    end

    """print(generate_eT_code_from_einsum(
        routine_name="$routine_name",
        prefactor=$prefactor,
        contraction_string="$einsum_str",
        arrays=[$tensor_str, $symbol],
        symbols=[$symbol_str, "$symbol"],
    ), end='!\\n!\\n')"""
end

# Print term in latex form
# TODO: Add possibility of changing tensors names
function print_latex(t::Term, translation, renaming)
    line = IOBuffer()
    # Sign
    if t.scalar > 0
        print(line, "+ ")
    else
        print(line, "- ")
    end
    # Prefactor
    if typeof(t.scalar) == Rational{Int64}
        if isone(denominator(t.scalar))
            if !isone(abs(numerator(t.scalar)))
                print(line, string(abs(numerator(t.scalar))) * " ")
            end
        else
            print(line, "\\frac{", string(abs(numerator(t.scalar))), "}{", string(denominator(t.scalar)), "} ")
        end
    elseif !isone(t.scalar)
        print(line, abs(string(t.scalar)) * " ")
    end

    translation = update_index_translation(t, translation)     # I have no idea if this line is needed or not. 
    # Sum
    if !isempty(t.sum_indices)
        print(line, "\\sum_{")
        for i in t.sum_indices
            print_mo_index(line, t.constraints, translation, i)
        end
        print(line, "} ")
    end
    # Deltas
    for d in t.deltas
        print(line, "\\delta_{")
        for p in d.indices
            print_mo_index(line, t.constraints, translation, p)
        end
        print(line, "} ")
    end
    # Tensors
    i = 1
    while i <= length(t.tensors)
        done = false
        j = 1
        while !done
            done = true
            if length(t.tensors) >= i + j && t.tensors[i+j] == t.tensors[i]
                j += 1
                done = false
            end
        end
        symb = get_symbol(t.tensors[i])
        for i in eachindex(renaming)
            if symb == renaming[i][1]
                symb = renaming[i][2]
            end
        end
        print(line, symb)
        inds = get_indices(t.tensors[i])
        if !isempty(inds)
            print(line, "_{")
            for ind in inds
                print_mo_index(line, t.constraints, translation, ind)
            end
            print(line, "}")
        end
        if j > 1
            print(line, "^{", string(j), "} ")
            i += j
        else
            print(line, " ")
            i += 1
        end
    end
    # Operators ##### LIMITED TO SINGLET EXCITATIONS OPERATORS
    for e in t.operators
        print(line, "E_{")
        print_mo_index(line, t.constraints, translation, e.p, e.q)
        print(line, "} ")
    end
    return String(take!(line))
end