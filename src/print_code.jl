export print_code
using Printf

print_code(t::Term) = print_code(t, "X")

function print_code(t::Term, symbol::String)
    # Print python np.einsum code
    scalar_str = @sprintf "%+12.8f" t.scalar

    external = prod(print_mo_index.(setdiff(get_all_indices(t), t.sum_indices)))

    # Make einsum_str
    einsum_str = "\""
    for a in t.tensors
        for b in get_indices(a)
            einsum_str *= print_mo_index(b)
        end
        einsum_str *= ","
    end
    einsum_str = einsum_str[begin:end-1] * "->" * external * "\""

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