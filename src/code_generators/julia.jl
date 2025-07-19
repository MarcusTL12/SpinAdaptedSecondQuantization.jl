export print_julia_function

"""
print_julia_function(name, ex::Expression, [symbol="X"];
                    [tensor_translation], [explicit_tensor_blocks])

Returns String with a runnable julia function to evaluate the given expression
`ex`. `name` specifies the name of the function, and the optional parameter
`symbol` specifies the name of the accumulator variable used inside the
function. If not specified it will use the generic name `X`.

The optional keyword argument `tensor_translation` is a dictionary for replacing
the name of certain tensors.

The optional keyword argument `explicit_tensor_blocks` is a vector of names of
tensors you want subblocks to be taken as separate inputs to the function
rather than using `@view` to get subblocks.

The returned function uses the
[`TensorOperations.jl`](https://github.com/Jutho/TensorOperations.jl)
package for tensor contractions.

# Example

```jldoctest
julia> using SpinAdaptedSecondQuantization

julia> L = psym_tensor("L", 1, 3, 4, 5);

julia> R = psym_tensor("R", 2, 3, 4, 5);

julia> ex = ∑(L * R * occupied(3, 5) * virtual(1, 2, 4), 3:5)
∑_icj(L_aicj R_bicj)
julia> print_julia_function("density", ex) |> print
function density(no, nv, L, R)
    nao = no + nv
    o = 1:no
    v = no+1:nao
    
    X = zeros(nv, nv)
    L_vovo = @view L[v,o,v,o]
    R_vovo = @view R[v,o,v,o]
    @tensoropt (a=>10χ,b=>10χ,c=>10χ,i=>χ,j=>χ) begin
        X[a,b] += L_vovo[a,i,c,j]*R_vovo[b,i,c,j]
    end
    X
end
julia> print_julia_function("density", ex; \
explicit_tensor_blocks=["L", "R"]) |> print
function density(no, nv, L_vovo, R_vovo)
    nao = no + nv
    o = 1:no
    v = no+1:nao
    
    X = zeros(nv, nv)
    @tensoropt (a=>10χ,b=>10χ,c=>10χ,i=>χ,j=>χ) begin
        X[a,b] += L_vovo[a,i,c,j]*R_vovo[b,i,c,j]
    end
    X
end
```

!!! note
    This function currently only supports expressions input and output tensors
    are indexed by occupied, virtual or general molecular orbitals.

    There is also not suppot for expressions containing Kronecker deltas, so
    these would require a bit of pre-processing.
"""
function print_julia_function(name, ex::Expression, symbol="X";
    tensor_translation=Dict(), explicit_tensor_blocks=String[])

    indices = get_external_indices(ex[1])
    external_constraints =
        Constraints([p => ex[1].constraints(p) for p in indices])

    for i in 2:length(ex.terms)
        if get_external_indices(ex[i]) != indices
            throw("Different external indices on different terms in same\
 expression is not supported for code generation!")
        end

        for p in indices
            if ex[i].constraints(p) != external_constraints(p)
                throw("Different index constraints on external indices on\
 different terms in same expression is not supported for code generation!")
            end
        end
    end

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

    out_dims = [space_dim_dict[external_constraints(p)] for p in indices]

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

    index_groups = Dict()

    for (p, s) in external_constraints
        push!(get!(index_groups, s, Int[]), p)
    end

    translation = translate(index_groups...)

    index_names = String[]
    io = IOBuffer()
    for p in indices
        print_mo_index(io, external_constraints, translation, p)
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

        local_trans = update_index_translation(t, IndexTranslation())

        term_indices = get_all_indices(t)

        for p in term_indices
            io = IOBuffer()
            print_mo_index(io, t.constraints, local_trans, p)
            n = String(take!(io))
            if !haskey(all_index_names, n)
                s = t.constraints(p)
                fac = if s ⊆ OccupiedOrbital
                    1
                elseif s ⊆ VirtualOrbital
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
