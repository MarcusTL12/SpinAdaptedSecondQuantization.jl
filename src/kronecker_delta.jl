export delta, δ

struct KroneckerDelta
    p::MOIndex
    q::MOIndex

    function KroneckerDelta(p::MOIndex, q::MOIndex)
        if isdisjoint(p, q)
            0
        elseif p == q
            1
        elseif p < q
            new(p, q)
        else
            new(q, p)
        end
    end
end

function Base.show(io::IO, d::KroneckerDelta)
    print(io, "δ_", d.p, d.q)
end

# Externally visible constructor
delta(p, q) = Expression(KroneckerDelta(p, q))

# Simple unicode alias
δ(p, q) = delta(p, q)

function Base.isless(d1::KroneckerDelta, d2::KroneckerDelta)
    (d1.p, d1.q) < (d2.p, d2.q)
end
