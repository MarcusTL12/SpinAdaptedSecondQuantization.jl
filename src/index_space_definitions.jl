export GeneralOrbital, OccupiedOrbital, VirtualOrbital

const GeneralOrbital = new_space(:GeneralOrbital, "g", "pqrstuv")
const VirtualOrbital = new_space(:VirtualOrbital, "v", "abcdefg")
const OccupiedOrbital = new_space(:OccupiedOrbital, "o", "ijklmno")

add_space_sum(OccupiedOrbital, VirtualOrbital, GeneralOrbital)

add_space_intersection(GeneralIndex, GeneralOrbital, GeneralOrbital)
add_space_intersection(GeneralOrbital, OccupiedOrbital, OccupiedOrbital)
add_space_intersection(GeneralOrbital, VirtualOrbital, VirtualOrbital)

# set_color(GeneralOrbital, :nothing)
# set_color(OccupiedOrbital, :light_green)
# set_color(VirtualOrbital, :cyan)
