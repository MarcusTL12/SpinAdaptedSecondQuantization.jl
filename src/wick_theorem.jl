# TODO T <: FermionOperator
function wick_theorem(opstring :: Vector{T}) where T <: Operator
    for op in opstring
        @assert typeof(op) == FermionOperator
    end
    list_of_pairs = fully_contracted_pairs(opstring)
    sum(prod(contract(a,b) for (a,b) in pairs) for pairs in list_of_pairs)
end

function contract(a :: FermionOperator, b :: FermionOperator)
    # Contractions used in Wick's theorem
    # contract(a, b) = <vac| a b |vac>
    δ(a.p, b.p) * (!a.dag && b.dag) * (a.spin == b.spin)
end

function fully_contracted_pairs(vec :: Vector{T}) where T
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
        x = fully_contracted_pairs(rest)
        for elem in x
            counter += 1
            list[counter] = vcat([pair], elem)
        end
    end

    return list
end