export print_eT_function_generator

"""
print_eT_function_generator(name, ex, symbol, indices, translation,
                    [wf_type], [tensor_translation], [noinput_tensors], [outperms])

"""
function print_eT_function_generator(name, ex::Expression, symbol, indices,
    translation, wf_type="ccs", tensor_translation=Dict(), noinput_tensors=String[],
    outperms=nothing)
    function get_tensor_name(t, n)
        get(tensor_translation, (t, n), get(tensor_translation, t, t))
    end

    function get_block_name(n, spaces)
        io = IOBuffer()
        print(io, "$(n)")
        if !isempty(spaces)
            print(io, "_")
            for s in spaces
                print(io, getshortname(s))
            end
        end
        String(take!(io))
    end

    function make_param_def(julianame, name, tens, spaces)
        perms = get_permutations(tens)

        filter!(perms) do perm
            permute!(copy(spaces), perm) == spaces
        end

        print_perms = length(perms) > 1 || !issorted(only(perms))

        is_inp = julianame ∉ noinput_tensors
        io = IOBuffer()
        if print_perms
            print(io, "$julianame = (\"$name\", $is_inp, $perms)")
        else
            print(io, "$julianame = (\"$name\", $is_inp)")
        end
        String(take!(io))
    end

    parameters = []

    outspaces = [translation(p)[1] for p in indices]

    outdims = [getshortname(s) for s in outspaces]

    func_body = IOBuffer()

    println(func_body, "func = FortranFunction((\"$symbol\", $outdims))")

    if !isnothing(outperms)
        println(func_body, "outperms = $outperms")
    end

    index_names = String[]
    io = IOBuffer()
    for p in indices
        print_mo_index(io, ex[1].constraints, translation, p)
        push!(index_names, String(take!(io)))
    end

    tensor_body = IOBuffer()

    for t in ex.terms
        if !isempty(t.operators)
            throw("Cannot make function to compute an operator!")
        end
        if !isempty(t.deltas)
            throw("Do not yet support generation of code with kronecker deltas!")
        end

        local_trans = update_index_translation(t, translation)

        print(tensor_body, "update_code!(func, ein\"")

        tensor_list = IOBuffer()

        print(tensor_list, "[")

        isfirst2 = true
        for tens in t.tensors
            if isfirst2
                isfirst2 = false
            else
                print(tensor_body, ",")
                print(tensor_list, ", ")
            end

            inds = get_indices(tens)

            spaces = [t.constraints(p) for p in inds]

            block_name = get_block_name(get_symbol(tens), spaces)
            tens_name = get_tensor_name(block_name, length(inds))
            julia_name = get_block_name(get_symbol(tens), spaces)

            param_def = make_param_def(julia_name, tens_name, tens, spaces)

            if param_def ∉ parameters
                push!(parameters, param_def)
            end

            print(tensor_list, julia_name)

            for p in inds
                print_mo_index(tensor_body, t.constraints, local_trans, p)
            end
        end

        print(tensor_body, "->")

        for p in index_names
            print(tensor_body, p)
        end

        print(tensor_list, "]")

        print(tensor_body, "\", ", t.scalar,
            ", ", String(take!(tensor_list)))

        if !isnothing(outperms)
            print(tensor_body, ", outperms")
        end

        println(tensor_body, ")")
    end

    for param in parameters
        println(func_body, param)
    end

    print(func_body, String(take!(tensor_body)))

    print(func_body, "finalize_eT_function(func, \"$name\", \"$wf_type\")")

    io = IOBuffer()

    println(io, "let")

    for l in eachsplit(String(take!(func_body)), "\n")
        println(io, "    ", l)
    end

    println(io, "end")
    String(take!(io))
end
