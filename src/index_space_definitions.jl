export GeneralOrbital, OccupiedOrbital, VirtualOrbital

const GeneralOrbital = new_space(:GeneralOrbital, "g", "pqrstuv")
const VirtualOrbital = new_space(:VirtualOrbital, "v", "abcdefg")
const OccupiedOrbital = new_space(:OccupiedOrbital, "o", "ijklmno")

add_space_sum(OccupiedOrbital, VirtualOrbital, GeneralOrbital)

add_subspace_relation(GeneralIndex, GeneralOrbital)
add_subspace_relation(GeneralOrbital, OccupiedOrbital)
add_subspace_relation(GeneralOrbital, VirtualOrbital)

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

# set_color(GeneralOrbital, :nothing)
# set_color(OccupiedOrbital, :light_green)
# set_color(VirtualOrbital, :cyan)
