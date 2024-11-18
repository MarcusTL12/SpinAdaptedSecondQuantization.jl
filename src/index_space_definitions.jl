export GeneralOrbital, OccupiedOrbital, VirtualOrbital

const GeneralOrbital = new_space(:GeneralOrbital, "g", "pqrstuv")
const VirtualOrbital = new_space(:VirtualOrbital, "v", "abcdefg")
const OccupiedOrbital = new_space(:OccupiedOrbital, "o", "ijklmno")

add_subspace_relation(GeneralIndex, GeneralOrbital)

add_space_sum(OccupiedOrbital, VirtualOrbital, GeneralOrbital)

export occupied, virtual, electron

function occupied(indices...)
    constrain(p => OccupiedOrbital for p in indices)
end

function virtual(indices...)
    constrain(p => VirtualOrbital for p in indices)
end

function electron(indices...)
    constrain(p => GeneralOrbital for p in indices)
end

const Boson = new_space(:Boson, "b", "αβγδ")
add_subspace_relation(GeneralIndex, Boson)

export boson_ind

function boson_ind(indices...)
    constrain(α => Boson for α in indices)
end
