
struct KroneckerDelta
    p :: MOIndex
    q :: MOIndex

    function KroneckerDelta(p::MOIndex, q::MOIndex)
        if isdisjoint(p, q)
            0
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
