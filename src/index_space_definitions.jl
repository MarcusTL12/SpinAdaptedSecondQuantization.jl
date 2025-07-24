export GeneralOrbital, OccupiedOrbital, VirtualOrbital

const GeneralOrbital = new_space(:GeneralOrbital, "g", "pqrstuv")
const VirtualOrbital = new_space(:VirtualOrbital, "v", "abcdefg")
const OccupiedOrbital = new_space(:OccupiedOrbital, "o", "ijklmno")

add_subspace_relation(GeneralIndex, GeneralOrbital)

add_space_sum(OccupiedOrbital, VirtualOrbital, GeneralOrbital)

export occupied, virtual, electron

"""
    occupied(indices...)

Constraining the given indices to the OccupiedOrbital index space.

# Example

```jldoctest
julia> using SpinAdaptedSecondQuantization

julia> occupied(1, 2)
C(i∈o, j∈o)

julia> E(1, 2) * ans
E_ij

```
"""
function occupied(indices...)
    constrain(p => OccupiedOrbital for p in indices)
end

"""
    virtual(indices...)

Constraining the given indices to the VirtualOrbital index space.

# Example

```jldoctest
julia> using SpinAdaptedSecondQuantization

julia> virtual(1, 2)
C(a∈v, b∈v)

julia> E(1, 2) * ans
E_ab

```
"""
function virtual(indices...)
    constrain(p => VirtualOrbital for p in indices)
end

"""
    electron(indices...)

Constraining the given indices to the GeneralOrbital index space.

# Example

```jldoctest
julia> using SpinAdaptedSecondQuantization

julia> electron(1, 2)
C(p∈g, q∈g)

julia> E(1, 2) * ans
E_pq

```
"""
function electron(indices...)
    constrain(p => GeneralOrbital for p in indices)
end
