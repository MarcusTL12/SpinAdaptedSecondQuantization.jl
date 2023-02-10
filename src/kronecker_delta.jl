export delta, δ

struct KroneckerDelta
    indices::Vector{MOIndex}

    function KroneckerDelta(indices)
        indices = unique!(sort(indices))

        if length(indices) <= 1
            1
        else
            for i in eachindex(indices), j in (i+1):length(indices)
                if isdisjoint(indices[i], indices[j])
                    return 0
                end
            end

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

function exchange_indices(d::KroneckerDelta, mapping)
    KroneckerDelta([exchange_index(p, mapping) for p in d.indices])
end
