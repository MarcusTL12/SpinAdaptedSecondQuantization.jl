function wick_theorem(string)
    for op in string
        @assert typeof(op) == SASQ.FermionOperator
    end

end

function contract(a, b)
    δ(a.p, b.p) * (!a.dag && b.dag) * (a.spin == b.spin)
end

function find_all_pairs(vec :: Vector{T}) where T
    # Recursive function which finds all possible unique pairs in a vector
    # [1,2] -> [[(1,2)]]
    # [1,2,3,4] ->  [[(1, 2), (3, 4)], [(1, 3), (2, 4)], [(1, 4), (2, 3)]]

    n = length(vec)
    if isodd(n)
        throw("length is odd")
    end
    if n == 2
        return [[(vec[1], vec[2])]]
    end

    # The number of ways to pair a list of n numbers is
    # the double factorial:
    # (n-1)!! = 1*3*5*...*(n-1).
    list = Vector{Vector{Tuple{T,T}}}(undef, prod(1:2:(n-1)))
    counter = 0
    for i = 2:n
        pair = (vec[1], vec[i])
        rest = vec[(1:n .!= 1) .&& (1:n .!= i)]
        x = find_all_pairs(rest)
        for elem in x
            counter += 1
            list[counter] = vcat([pair], elem)
        end
    end

    return list
end

for i = 2:2:16
    v = collect(1:i)
    @show length(find_all_pairs(v))
end