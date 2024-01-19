export print_eT_function_generator

"""
print_eT_function_generator(name, ex, symbol, indices, translation,
                    [wf_type], [tensor_translation], [noinput_tensors])

"""
function print_eT_function_generator(name, ex::Expression, symbol, indices,
    translation, wf_type="ccs", tensor_translation=Dict(), noinput_tensors=String[])
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

    space_dim_dict = Dict([
        OccupiedOrbital => "o",
        VirtualOrbital => "v",
        GeneralOrbital => "g",
    ])

    function make_param_def(tens_symb, julianame, name, spaces, flag)
        io = IOBuffer()
        if isempty(spaces)
            funcname = if flag == :I && tens_symb ∉ noinput_tensors
                "input_scalar"
            elseif flag == :O
                "output_scalar"
            else
                "Sym"
            end
            print(io, "$julianame = $funcname(\"$name\")")
        else
            funcname = if flag == :I && tens_symb ∉ noinput_tensors
                "input_tensor"
            elseif flag == :O
                "output_tensor"
            else
                "noio_tensor"
            end
            print(io, "$julianame = $funcname(\"$name\" => (")
            isfirst = true
            for space in spaces
                if isfirst
                    isfirst = false
                else
                    print(io, ", ")
                end
                print(io, '"', space_dim_dict[space], '"')
            end

            print(io, "))")
        end

        String(take!(io))
    end

    parameters = []

    outspaces = [translation(p)[1] for p in indices]

    func_body = IOBuffer()

    index_names = String[]
    io = IOBuffer()
    for p in indices
        print_mo_index(io, ex[1].constraints, translation, p)
        push!(index_names, String(take!(io)))
    end

    all_index_names = Dict()

    tensor_body = IOBuffer()

    for t in ex.terms
        if !isempty(t.operators)
            throw("Cannot make function to compute an operator!")
        end
        if !isempty(t.deltas)
            throw("Do not yet support generation of code with kronecker deltas!")
        end

        local_trans = update_index_translation(t, translation)

        term_indices = get_all_indices(t)

        for p in term_indices
            io = IOBuffer()
            print_mo_index(io, t.constraints, local_trans, p)
            n = String(take!(io))
            if !haskey(all_index_names, n)
                s = t.constraints(p)
                fac = if s <: OccupiedOrbital
                    1
                elseif s <: VirtualOrbital
                    10
                else
                    11
                end

                all_index_names[n] = fac
            end
        end

        print(tensor_body, "    $symbol")
        isfirst = true
        print(tensor_body, "[")
        for x in index_names
            if isfirst
                isfirst = false
            else
                print(tensor_body, ",")
            end
            print(tensor_body, "$x")
        end
        print(tensor_body, "]")

        print(tensor_body, " += ")

        isfirst2 = true

        if isone(-t.scalar)
            print(tensor_body, "-")
        elseif !isone(t.scalar)
            printscalar(tensor_body, t.scalar)
            isfirst2 = false
        end

        for tens in t.tensors
            if isfirst2
                isfirst2 = false
            else
                print(tensor_body, "*")
            end

            inds = get_indices(tens)

            spaces = [t.constraints(p) for p in inds]

            tens_name = get_tensor_name(get_symbol(tens), length(inds))
            julia_name = get_block_name(get_symbol(tens), spaces)
            block_name = get_block_name(tens_name, spaces)

            param_def = make_param_def(get_symbol(tens), julia_name, block_name, spaces, :I)

            if param_def ∉ parameters
                push!(parameters, param_def)
            end

            print(tensor_body, julia_name)

            isfirstinner = true
            for p in inds
                if isfirstinner
                    print(tensor_body, "[")
                    isfirstinner = false
                else
                    print(tensor_body, ",")
                end
                print_mo_index(tensor_body, t.constraints, local_trans, p)
            end
            if !isfirstinner
                print(tensor_body, "]")
            end
        end

        println(tensor_body, "")
    end

    println(func_body, "reset_state()")

    println(func_body, make_param_def(symbol, symbol, symbol, outspaces, :O))

    for param in parameters
        println(func_body, param)
    end

    all_index_names = sort!(collect(all_index_names))

    print(func_body, "@tensor opt=(")
    isfirst = true
    for (n, f) in all_index_names
        if isfirst
            isfirst = false
        else
            print(func_body, ",")
        end
        print(func_body, "$n=>")
        if f != 1
            print(func_body, f)
        end
        print(func_body, "χ")
    end
    println(func_body, ") backend=eTbackend begin")

    print(func_body, String(take!(tensor_body)))
    println(func_body, "end")

    print(func_body, "finalize_eT_function(\"$name\", \"$wf_type\")")

    io = IOBuffer()

    println(io, "let")

    for l in eachsplit(String(take!(func_body)), "\n")
        println(io, "    ", l)
    end

    println(io, "end")
    String(take!(io))
end
