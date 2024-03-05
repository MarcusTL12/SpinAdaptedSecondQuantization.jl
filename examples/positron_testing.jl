include("multilevel.jl")

right_state = fermiondag(1, α) * ivir(1)

function act_on_actual_ket(O)
    act_on_ket(O * right_state)
end


using Printf
function print_code_einsum_testing(t::SASQ.Term, symbol::String, translation, fixed)
    # Print python np.einsum code for term t added to symbol. With deltas, only works on 0,2 or 4 dimensional symbol
    # fixed = [] or fixed = ["a", "i"]
    scalar_str = @sprintf "%+12.8f" t.scalar
    translation = SASQ.update_index_translation(t, translation)

    function write_extract(t, a, external, translation)
        # t term, a tensor
        write_str = " extract_mat($(SASQ.get_symbol(a)), \""
        for b in SASQ.get_indices(a)
            if b ∈ t.sum_indices || sprint(SASQ.print_mo_index, t.constraints, translation, b)[1] in external
                if  t.constraints[b] ⊆ VirtualOrbital
                    write_str *= t.constraints[b] ⊆ ActiveOrbital ? "v" : "V"
                else
                    write_str *= "o"
                end
            else
                #Here only if is a positron not summed and not a real external
                write_str *= b == 1 ? "I" : "A"
            end
        end
        return write_str * "\", o, v)"
    end

    # Finding deltas. Only δ_ab and/or δ_ac. No δ_bc
    fix_b = false
    fix_c = false
    fix_j = false
    fix_k = false
    if length(t.deltas) > 0
        for d in t.deltas
            delta_ind = sprint(SASQ.print_mo_index, t.constraints, translation, d.indices...)
            if 'b' in delta_ind
                fix_b = true
            end
            if 'c' in delta_ind
                fix_c = true
            end
            if 'j' in delta_ind
                fix_j = true
            end
            if 'k' in delta_ind
                fix_k = true
            end
        end
    end

    external_int = SASQ.get_external_indices(t)
    external = sprint(SASQ.print_mo_index, t.constraints, translation, external_int...)

    # Remove a and i  from external
    external = join(reverse([a for a in external if a ∉ fixed]))
    new_ext = external

    # Determining actual externals, after the deltas

    # No deltas in this case
    # new_ext = ""
    # if length(external) >= 2
    #     new_ext *= fix_b ? "" : "b"
    #     new_ext *= fix_j ? "" : "j"
    # end
    # if length(external) == 4
    #     new_ext *= fix_c ? "" : "c"
    #     new_ext *= fix_k ? "" : "k"
    # end

    # temp_color = index_color
    # disable_color()

    # Make einsum_str and find not-summed tensors
    einsum_str = "\""
    notsum_str = ""
    not_summed_tensors = []
    print_einsum = false
    for a in t.tensors
        # if (No indices) or (nothing summed  and 0 external)                  
        if length(SASQ.get_indices(a)) == 0 || (length(t.sum_indices) == 0 && sum([1 for b in SASQ.get_indices(a) if sprint(SASQ.print_mo_index, t.constraints, translation, b)[1] in new_ext]) == 0)
            push!(not_summed_tensors, a)
        else
            indices = []
            for b in SASQ.get_indices(a)
                if b ∈ t.sum_indices || sprint(SASQ.print_mo_index, t.constraints, translation, b)[1] in new_ext
                    push!(indices, sprint(SASQ.print_mo_index, t.constraints, translation, b))
                end
            end
            if join(indices) == new_ext
                push!(not_summed_tensors, a)
                new_ext = ""
            elseif join(indices) == ""
                    push!(not_summed_tensors, a)
            else
                einsum_str *= join(indices)
                einsum_str *= ","
                print_einsum = true
            end
        end
    end

    if print_einsum
        einsum_str = einsum_str[begin:end-1] * "->" * new_ext * "\""
    end

    # if (temp_color)
    #     enable_color()
    # end

    # Make tensor_str
    tensor_str = ""
    for a in t.tensors 
        if a ∉ not_summed_tensors
            tensor_str *= ","
            tensor_str *= write_extract(t, a, external, translation)
        end
    end

    # Make string for not-summed tensors
    for a in not_summed_tensors
        if length(SASQ.get_indices(a)) > 0
            notsum_str *= " *"
            notsum_str *= write_extract(t, a, external, translation)
        else
            notsum_str *= " * $(SASQ.get_symbol(a))"
        end
    end

    # Printing depending on deltas
    if length(external) >= 2
        pre_string = "$(symbol)_$(external)"
        pre_string *= fix_b ? "[a" : "[:"
        pre_string *= fix_j ? ",i" : ",:"
        if length(external) == 4
            pre_string *= fix_c ? ",a" : ",:"
            pre_string *= fix_k ? ",i" : ",:"
        end
        pre_string *= "]"
    else
        pre_string = "$(symbol)"
    end

    if length(tensor_str) > 0
        return "$pre_string = $pre_string .+ $scalar_str $notsum_str * np.einsum($einsum_str$tensor_str, optimize=\"optimal\");"
    else
        return "$pre_string += $scalar_str $notsum_str;"
    end
end

function contains_delta(term)
    # println(SASQ.get_external_indices(term))
    if length(term.deltas) == 0
        return true
    end
    return false
end

function contains_wrong_s(term)
    # Check if there is s_ai in that term
    for tens in term.tensors
        if ((tens.indices[1] == 1 || tens.indices[2] == 2) && tens.symbol in ["s","s2"])
            return false
        end
    end
    return true
end

function filter_unwanted(expression)
    terms = [contains_delta(t) && contains_wrong_s(t) for t in expression.terms]
    return SASQ.Expression(expression[terms])
end

h_elec = ∑(real_tensor("h", 1, 2) * E(1, 2) * active(1, 2), 1:2)
h_posi = ∑(real_tensor("h_p", 1, 2) * E(1, 2) * ivir(1, 2), 1:2)

h = h_elec + h_posi

g_elec = 1 // 2 * ∑(psym_tensor("g", 1, 2, 3, 4) * e(1, 2, 3, 4) * active(1, 2, 3, 4), 1:4) |> simplify
g_int = -∑(psym_tensor("g_p", 1, 2, 3, 4) * e(1, 2, 3, 4) * active(1, 2) * ivir(3, 4), 1:4) |> simplify
# g_posi = 1 // 2 * ∑(psym_tensor("g_pp", 1, 2, 3, 4) * e(1, 2, 3, 4) * ivir(1, 2, 3, 4), 1:4) |> simplify

H_elec = h_elec + g_elec
H = H_elec + h_posi + g_int

# E_hf = hf_expectation_value(right_state' * H * right_state) |> simplify_heavy
# @show E_hf
# println()

hF_elec = ∑((real_tensor("F", 1, 2) +
             ∑(aocc(3) * (-2psym_tensor("g", 1, 2, 3, 3) +
                          psym_tensor("g", 1, 3, 3, 2)), [3])) * E(1, 2) * active(1, 2), 1:2)

HF_elec = simplify(hF_elec + g_elec)

HF = HF_elec + h_posi + g_int

E_hf = hf_expectation_value(right_state' * (H + HF) // 2 * right_state) |> simplify_heavy
E_hf = filter_unwanted(E_hf)

@show E_hf  # Filtering is not needed

open("file_energy.py", "w") do output_file
    for t in E_hf.terms
        println(output_file, print_code_einsum_testing(t, "E", SASQ.IndexTranslation(), ['I']))
    end
end

println()

ex_ketop(a, i) = E(a, i) * aocc(i) * avir(a)
ex_ketop(a, i, b, j) = E(a, i) * E(b, j) * aocc(i, j) * avir(a, b)

deex_braop(a, i) = 1 // 2 * ex_ketop(a, i)'
deex_braop(a, i, b, j) = 1 // 3 * ex_ketop(a, i, b, j)' +
                         1 // 6 * ex_ketop(a, j, b, i)'

ex_positron(a, b) = E(a, b) * ivir(a, b)
ex_positron(a, b, c, d) = E(a, i) * E(b, j) * ivir(a, b, c, d)

deex_positron(a, i) = 1 // 2 * ex_positron(a, i)'     #probabily this 1/2 can be removed
deex_positron(a, i, b, j) = 1 // 3 * ex_positron(a, i, b, j)' +
                            1 // 6 * ex_positron(a, j, b, i)'

t(inds...) = psym_tensor("t", inds...)
s(inds...) = psym_tensor("s", inds...)

T2 = 1 // 2 * ∑(
    t(1:4...) * ex_ketop(1, 2, 3, 4),
    1:4
)

S1 = ∑(s(1, 2, 3, 4) * ex_ketop(1, 2) * ex_positron(3, 4), 1:4)
S2 = 1 // 2 * ∑(
    psym_tensor("s2", 1:6...) * ex_ketop(1, 2, 3, 4) * ex_positron(5, 6),
    1:6
)

T = T2 + S1 + S2

@show HF

function omega(proj, op, n)
    hf_expectation_value(simplify(right_state' * proj * bch(op, T, n) * right_state))
end

function omega_AI()
    o = omega(deex_positron(2, 1), HF, 1)
    o = simplify_heavy(o)
    o = look_for_tensor_replacements_smart(o, make_exchange_transformer("t", "u"))
    o = look_for_tensor_replacements_smart(o, make_exchange_transformer("g", "L"))
    return filter_unwanted(o)
end

Omega_AI = omega_AI()

open("file_omega_AI.py", "w") do output_file
    for t in Omega_AI.terms
        println(output_file, print_code_einsum_testing(t, "Omega_AI", SASQ.IndexTranslation(), ['A','I']))
    end
end

function omega_ai()
    o = omega(deex_braop(4,3), HF, 1)
    o = simplify_heavy(o)
    o = look_for_tensor_replacements_smart(o, make_exchange_transformer("t", "u"))
    o = look_for_tensor_replacements_smart(o, make_exchange_transformer("g", "L"))
    return filter_unwanted(o)
end

# ...

@show omega_AI()
@show omega_ai()
;
