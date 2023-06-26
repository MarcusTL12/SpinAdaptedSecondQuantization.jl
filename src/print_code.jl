export print_code
using Printf

""" Print term to numpy's einsum routine
"""
print_code(t::Term, translation) = print_code(t, "X", translation)

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

    println("$(symbol)_$(external) += $scalar_str * np.einsum($einsum_str$tensor_str, optimize=\"optimal\");")
end

""" Used with https://github.com/alexancp/einsumpath-to-eT to generate Fortran subroutines for eT.
"""
function print_eT_code(t :: Term, symbol, translation, routine_name)
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
