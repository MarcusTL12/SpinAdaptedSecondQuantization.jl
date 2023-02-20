export delta, δ

struct KroneckerDelta
    indices::Vector{MOIndex}

    function KroneckerDelta(indices)
        indices = unique!(sort(indices))

        if length(indices) <= 1
            1
        else
            new(indices)
        end
    end
end

function KroneckerDelta(indices...)
    KroneckerDelta(collect(indices))
end

function Base.show(io::IO, d::KroneckerDelta)
    print(io, "δ_")

    for p in d.indices
        print(io, p)
    end
end

# Externally visible constructor
delta(indices...) = Expression(KroneckerDelta(indices...))

# Simple unicode alias
δ(indices...) = delta(indices...)

function Base.isless(d1::KroneckerDelta, d2::KroneckerDelta)
    d1.indices < d2.indices
end

function Base.:(==)(d1::KroneckerDelta, d2::KroneckerDelta)
    d1.indices == d2.indices
end

function exchange_indices(d::KroneckerDelta, mapping)
    KroneckerDelta([exchange_index(p, mapping) for p in d.indices])
end

function compact_deltas(deltas::Vector{KroneckerDelta})
    forest = DisjointSets{MOIndex}()

    function add_index(i)
        if !haskey(forest.intmap, i)
            push!(forest, i)
        end
    end

    function add_index(r, i)
        if !haskey(forest.intmap, i)
            push!(forest, i)
        end

        union!(forest, r, i)
    end

    for d in deltas
        r, rest = Iterators.peel(d.indices)
        add_index(r)

        for i in rest
            add_index(r, i)
        end
    end

    group_dict = Dict{MOIndex,Vector{MOIndex}}()

    for i in collect(forest)
        r = find_root!(forest, i)
        if haskey(group_dict, r)
            push!(group_dict[r], i)
        else
            group_dict[r] = [i]
        end
    end

    new_deltas = KroneckerDelta[]

    for v in values(group_dict)
        d = KroneckerDelta(v)
        if d == 0
            return 0
        elseif d == 1
            @warn "Fusing deltas should not be able to produce 1!"
        else
            push!(new_deltas, d)
        end
    end

    sort!(new_deltas)
end

#function space(d::KroneckerDelta)
#    s, rest = Iterators.peel(space(p) for p in d.indices)
#
#    for s2 in rest
#        s = typeintersect(s, s2)
#    end
#
#    s
#end
#