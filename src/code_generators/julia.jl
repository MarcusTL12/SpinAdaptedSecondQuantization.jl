export print_julia_function

"""
print_julia_function(name, ex, symbol, indices, translation,
                    [tensor_translation], [explicit_tensor_blocks])

"""
function print_julia_function(name, ex::Expression, symbol, indices, translation,
    tensor_translation=Dict(), explicit_tensor_blocks=String[])
    function get_tensor_name(t, n)
        get(tensor_translation, (t, n), get(tensor_translation, t, t))
    end

    function get_block_name(n, spaces)
        io = IOBuffer()
        print(io, "$(n)_")
        for s in spaces
            print(io, getshortname(s))
        end
        String(take!(io))
    end

    parameters = ["no", "nv"]

    space_dim_dict = Dict([
        OccupiedOrbital => "no",
        VirtualOrbital => "nv",
        GeneralOrbital => "nao",
    ])

    space_range_dict = Dict([
        OccupiedOrbital => "o",
        VirtualOrbital => "v",
        GeneralOrbital => ":",
    ])

    out_dims = [space_dim_dict[translation(p)[1]] for p in indices]

    func_body = IOBuffer()

    print(
        func_body,
        """
nao = no + nv
o = 1:no
v = no+1:nao

$symbol = """)

    if isempty(out_dims)
        println(func_body, "0.0")
    else
        isfirst = true
        for d in out_dims
            if isfirst
                isfirst = false
                print(func_body, "zeros(")
            else
                print(func_body, ", ")
            end
            print(func_body, "$d")
        end
        println(func_body, ")")
    end

    tensor_blocks = Dict()

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
        if !isempty(index_names)
            isfirst = true
            for x in index_names
                if isfirst
                    isfirst = false
                    print(tensor_body, "[")
                else
                    print(tensor_body, ",")
                end
                print(tensor_body, "$x")
            end
            print(tensor_body, "]")
        end

        print(tensor_body, " += ")

        isfirst = true

        if isone(-t.scalar)
            print(tensor_body, "-")
        elseif !isone(t.scalar)
            printscalar(tensor_body, t.scalar)
            isfirst = false
        end

        for tens in t.tensors
            if isfirst
                isfirst = false
            else
                print(tensor_body, "*")
            end

            inds = get_indices(tens)

            spaces = [t.constraints(p) for p in inds]

            tens_name = get_tensor_name(get_symbol(tens), length(inds))
            block_name = get_block_name(tens_name, spaces)

            if tens_name ∈ explicit_tensor_blocks
                if block_name ∉ parameters
                    push!(parameters, block_name)
                end
            else
                if tens_name ∉ parameters
                    push!(parameters, tens_name)
                end

                blocks = get(tensor_blocks, tens_name, [])
                if spaces ∉ blocks
                    push!(blocks, spaces)
                end
                tensor_blocks[tens_name] = blocks
            end

            print(tensor_body, block_name)

            isfirstinner = true
            for p in inds
                if isfirstinner
                    isfirstinner = false
                    print(tensor_body, "[")
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

    for (tens_name, blocks) in tensor_blocks
        for spaces in blocks
            print(func_body, get_block_name(tens_name, spaces))
            print(func_body, " = @view ", tens_name)
            isfirst = true
            for space in spaces
                if isfirst
                    isfirst = false
                    print(func_body, "[")
                else
                    print(func_body, ",")
                end
                print(func_body, space_range_dict[space])
            end
            println(func_body, "]")
        end
    end

    all_index_names = sort!(collect(all_index_names))

    print(func_body, "@tensoropt (")
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
    println(func_body, ") begin")

    print(func_body, String(take!(tensor_body)))
    print(func_body, "end")

    io = IOBuffer()

    print(io, "function ", name)
    isfirst = true
    for param in parameters
        if isfirst
            isfirst = false
            print(io, "(")
        else
            print(io, ", ")
        end
        print(io, param)
    end
    println(io, ")")

    for l in eachsplit(String(take!(func_body)), "\n")
        println(io, "    ", l)
    end

    println(io, "    $symbol\nend")
    String(take!(io))
end
