export print_code, print_eT_code
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

function print_code_einsum_withextract(t::Term, symbol::String, translation)
    # Print python np.einsum code
    scalar_str = @sprintf "%+12.8f" t.scalar
    translation = update_index_translation(t, translation)

    fix_a = false
    fix_i = false
    if length(t.deltas) > 0
        for d in t.deltas
            if 'a' in sprint(SASQ.print_mo_index, t.constraints, translation, d.indices...)
                fix_a = true
            elseif 'i' in sprint(SASQ.print_mo_index, t.constraints, translation, d.indices...)
                fix_i = true
            else
                throw("Unsimplified deltas")
            end
        end
    end

    external_int = get_external_indices(t)
    external = sprint(SASQ.print_mo_index, t.constraints, translation, external_int...)

    # Remove a and i  from external
    fixed = ['a', 'i']
    external = join([a for a in external if a ∉ fixed])

    temp_color = index_color
    disable_color()

    # Make einsum_str and find not-summed tensors
    einsum_str = "\""
    notsum_str = ""
    not_summed_tensors = []
    print_einsum = false
    for a in t.tensors
        # indices = [ind for ind in get_indices(a) if ind in t.sum_indices]
        indices = []
        for b in get_indices(a)
            if b ∈ t.sum_indices || sprint(SASQ.print_mo_index, t.constraints, translation, b)[1] in external
                push!(indices, sprint(SASQ.print_mo_index, t.constraints, translation, b))
            end
        end

        if length(get_indices(a)) > 0 && (length(t.sum_indices) > 0 || sum([1 for b in get_indices(a) if sprint(SASQ.print_mo_index, t.constraints, translation, b)[1] in external]) == 1)
            einsum_str *= join(indices)
            einsum_str *= ","
            print_einsum = true
        else
            push!(not_summed_tensors, a)
        end
    end
    
    new_ext = external
    if fix_a
        new_ext = "j"
        if fix_i
            new_ext = ""
        end
    elseif fix_i
        new_ext = "b"
    end

    if print_einsum
        einsum_str = einsum_str[begin:end-1] * "->" * new_ext * "\""
    end

    if (temp_color)
        enable_color()
    end

    # Make tensor_str
    tensor_str = ""
    for a in t.tensors 
        if a ∉ not_summed_tensors
            tensor_str *= ", extract_mat($(get_symbol(a)), \""
            for b in get_indices(a)
                if t.constraints[b] == VirtualOrbital
                    if b ∈ t.sum_indices || sprint(SASQ.print_mo_index, t.constraints, translation, b)[1] in external
                        tensor_str *= "v"
                    else
                        tensor_str *= "a"
                    end
                elseif t.constraints[b] == OccupiedOrbital
                    if b ∈ t.sum_indices || sprint(SASQ.print_mo_index, t.constraints, translation, b)[1] in external
                        tensor_str *= "o"
                    else
                        tensor_str *= "i"
                    end
                else
                    throw("Space not supported")
                end
            end
            tensor_str *= "\", o, v)"
        end
    end

    # Make string for not-summed tensors
    for a in not_summed_tensors
        if length(get_indices(a)) > 0
            notsum_str *= " * extract_mat($(get_symbol(a)), \""
            for b in get_indices(a)
                if t.constraints[b] == VirtualOrbital
                    if sprint(SASQ.print_mo_index, t.constraints, translation, b)[1] in external
                        notsum_str *= "v"
                    else
                        notsum_str *= "a"
                    end
                elseif t.constraints[b] == OccupiedOrbital
                    if sprint(SASQ.print_mo_index, t.constraints, translation, b)[1] in external
                        notsum_str *= "o"
                    else
                        notsum_str *= "i"
                    end
                else
                    throw("Space not supported")
                end
            end
            notsum_str *= "\", o, v)"
        else
            notsum_str *= " * $(get_symbol(a))"
        end
    end

    # Printing depending on sum or not
    pre_string = "$(symbol)_$(external)"
    if fix_a
        pre_string *= "[a,"
        if fix_i
            pre_string *= "i]"
        else
            pre_string *= ":]"
        end
    elseif fix_i
        pre_string *= "[:,i]"
    end

    if length(external) == 0
        pre_string = "$(symbol)"
    end

    if length(tensor_str) > 0
        return "$pre_string = $pre_string .+ $scalar_str $notsum_str * np.einsum($einsum_str$tensor_str, optimize=\"optimal\");"
    else
        return "$pre_string += $scalar_str $notsum_str;"
    end

end

function print_code_einsum_withextract_notfixed(t::Term, symbol::String, translation)
    # Print python np.einsum code
    scalar_str = @sprintf "%+12.8f" t.scalar
    translation = update_index_translation(t, translation)

    fix_a = false
    fix_i = false
    if length(t.deltas) > 0
        for d in t.deltas
            if 'a' in sprint(SASQ.print_mo_index, t.constraints, translation, d.indices...)
                fix_a = true
            elseif 'i' in sprint(SASQ.print_mo_index, t.constraints, translation, d.indices...)
                fix_i = true
            else
                throw("Unsimplified deltas")
            end
        end
    end

    external_int = get_external_indices(t)
    external = sprint(SASQ.print_mo_index, t.constraints, translation, external_int...)

    # Remove a and i  from external
    fixed = []
    external = join([a for a in external if a ∉ fixed])

    temp_color = index_color
    disable_color()

    # Make einsum_str and find not-summed tensors
    einsum_str = "\""
    notsum_str = ""
    not_summed_tensors = []
    print_einsum = false
    for a in t.tensors
        # indices = [ind for ind in get_indices(a) if ind in t.sum_indices]
        indices = []
        for b in get_indices(a)
            if b ∈ t.sum_indices || sprint(SASQ.print_mo_index, t.constraints, translation, b)[1] in external
                push!(indices, sprint(SASQ.print_mo_index, t.constraints, translation, b))
            end
        end

        if length(get_indices(a)) > 0 && (length(t.sum_indices) > 0 || sum([1 for b in get_indices(a) if sprint(SASQ.print_mo_index, t.constraints, translation, b)[1] in external]) == 1)
            einsum_str *= join(indices)
            einsum_str *= ","
            print_einsum = true
        else
            push!(not_summed_tensors, a)
        end
    end
    
    new_ext = external
    if fix_a
        new_ext = "j"
        if fix_i
            new_ext = ""
        end
    elseif fix_i
        new_ext = "b"
    end

    if print_einsum
        einsum_str = einsum_str[begin:end-1] * "->" * new_ext * "\""
    end

    if (temp_color)
        enable_color()
    end

    # Make tensor_str
    tensor_str = ""
    for a in t.tensors 
        if a ∉ not_summed_tensors
            tensor_str *= ", extract_mat($(get_symbol(a)), \""
            for b in get_indices(a)
                if t.constraints[b] == VirtualOrbital
                    if b ∈ t.sum_indices || sprint(SASQ.print_mo_index, t.constraints, translation, b)[1] in external
                        tensor_str *= "v"
                    else
                        tensor_str *= "a"
                    end
                elseif t.constraints[b] == OccupiedOrbital
                    if b ∈ t.sum_indices || sprint(SASQ.print_mo_index, t.constraints, translation, b)[1] in external
                        tensor_str *= "o"
                    else
                        tensor_str *= "i"
                    end
                else
                    throw("Space not supported")
                end
            end
            tensor_str *= "\", o, v)"
        end
    end

    # Make string for not-summed tensors
    for a in not_summed_tensors
        if length(get_indices(a)) > 0
            notsum_str *= " * extract_mat($(get_symbol(a)), \""
            for b in get_indices(a)
                if t.constraints[b] == VirtualOrbital
                    if sprint(SASQ.print_mo_index, t.constraints, translation, b)[1] in external
                        notsum_str *= "v"
                    else
                        notsum_str *= "a"
                    end
                elseif t.constraints[b] == OccupiedOrbital
                    if sprint(SASQ.print_mo_index, t.constraints, translation, b)[1] in external
                        notsum_str *= "o"
                    else
                        notsum_str *= "i"
                    end
                else
                    throw("Space not supported")
                end
            end
            notsum_str *= "\", o, v)"
        else
            notsum_str *= " * $(get_symbol(a))"
        end
    end

    # Printing depending on sum or not
    pre_string = "$(symbol)_$(external)"
    if fix_a
        pre_string *= "[a,"
        if fix_i
            pre_string *= "i]"
        else
            pre_string *= ":]"
        end
    elseif fix_i
        pre_string *= "[:,i]"
    end

    if length(external) == 0
        pre_string = "$(symbol)"
    end

    if length(tensor_str) > 0
        return "$pre_string = $pre_string .+ $scalar_str $notsum_str * np.einsum($einsum_str$tensor_str, optimize=\"optimal\");"
    else
        return "$pre_string += $scalar_str $notsum_str;"
    end

end