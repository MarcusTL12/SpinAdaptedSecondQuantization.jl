export wick_theorem

# TODO T <: FermionOperator
function wick_theorem(ex :: Expression)
    # Returns an Expression for <ref| ex1 |ref>
    sum(wick_theorem(term) for term in ex.terms)
end

function wick_theorem(t :: Term)
    # Returns an Expression for <ref| term |ref>
    wick_expr =  wick_theorem(t.operators)
    noop_term = noop_part(t)
    Expression([fuse(noop_term, wterm) for wterm in wick_expr.terms])
end

function wick_theorem(opstring :: Vector{T}) where T <: Operator
    # Returns an Expression for <ref| opstring |ref>
    for op in opstring
        @assert typeof(op) == FermionOperator
    end

    if isodd(length(opstring)) || length(opstring) == 0
        return Expression(0)
    end

    list_of_pairs = fully_contracted_pairs(opstring)
    list_of_pairs_index = fully_contracted_pairs(collect(1:length(opstring)))
    signs = map(find_sign, list_of_pairs_index)

    sum(signs[i] * prod(contract(a,b) for (a,b) in pairs) for (i, pairs) in enumerate(list_of_pairs))
end

function contract(a :: FermionOperator, b :: FermionOperator)
    # Contractions used in Wick's theorem
    # contract(a, b) = <ref| a b |ref>
    δ(a.p, b.p) * (a.dag && !b.dag) * (a.spin == b.spin) * occupied(a.p) +
    δ(a.p, b.p) * (!a.dag && b.dag) * (a.spin == b.spin) * virtual(a.p)
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

function find_sign(pairs)
    # Wick's theorem to find sign
    # sign = (-1)^number of intersections between pairs
    C = 1
    for i = 1:length(pairs), j = i+1:length(pairs)
        p1 = pairs[i]
        p2 = pairs[j]

        inter = intersect(p1[1]:p1[2], p2[1]:p2[2])
        if !isempty(inter) && inter != p1[1]:p1[2] && inter != p2[1]:p2[2]
            C *= -1
        end
    end

    return C
end