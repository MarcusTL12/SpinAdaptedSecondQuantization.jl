export GeneralOrbital, OccupiedOrbital, VirtualOrbital

const GeneralOrbital = :GeneralOrbital
const OccupiedOrbital = :OccupiedOrbital
const VirtualOrbital = :VirtualOrbital

add_space_names(GeneralOrbital, "g", "pqrstuv")
add_space_names(OccupiedOrbital, "o", "ijklmno")
add_space_names(VirtualOrbital, "v", "abcdefg")

add_space_sum(OccupiedOrbital, VirtualOrbital, GeneralOrbital)

add_space_intersection(GeneralIndex, GeneralOrbital, GeneralOrbital)
add_space_intersection(GeneralOrbital, OccupiedOrbital, OccupiedOrbital)
add_space_intersection(GeneralOrbital, VirtualOrbital, VirtualOrbital)

# set_color(GeneralOrbital, :nothing)
# set_color(OccupiedOrbital, :light_green)
# set_color(VirtualOrbital, :cyan)
