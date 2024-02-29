using SpinAdaptedSecondQuantization

# Define spaces with names and colors

ActiveOrbital = SASQ.new_space(:ActiveOrbital, "a", "pqrstuv")
InactiveOrbital = SASQ.new_space(:InactiveOrbital, "i", "pqrstuv")

set_color(ActiveOrbital, :cyan)
set_color(InactiveOrbital, :red)

ActiveVirtualOrbital = SASQ.new_space(:ActiveVirtualOrbital, "av", "abcdefg")
ActiveOccupiedOrbital = SASQ.new_space(:ActiveOccupiedOrbital, "ao", "ijklmno")

set_color(ActiveVirtualOrbital, :cyan)
set_color(ActiveOccupiedOrbital, :cyan)

InactiveVirtualOrbital = SASQ.new_space(:InactiveVirtualOrbital, "iv", "abcdefg")
InactiveOccupiedOrbital = SASQ.new_space(:InactiveOccupiedOrbital, "io", "ijklmno")

set_color(InactiveVirtualOrbital, :red)
set_color(InactiveOccupiedOrbital, :red)

# Define relations

SASQ.add_space_sum(ActiveOrbital, InactiveOrbital, GeneralOrbital)

SASQ.add_subspace_relation(GeneralOrbital, ActiveOrbital)
SASQ.add_subspace_relation(GeneralOrbital, InactiveOrbital)


SASQ.add_space_sum(ActiveVirtualOrbital, InactiveVirtualOrbital, VirtualOrbital)

SASQ.add_subspace_relation(VirtualOrbital, ActiveVirtualOrbital)
SASQ.add_subspace_relation(VirtualOrbital, InactiveVirtualOrbital)


SASQ.add_space_sum(ActiveOccupiedOrbital, InactiveOccupiedOrbital, OccupiedOrbital)

SASQ.add_subspace_relation(OccupiedOrbital, ActiveOccupiedOrbital)
SASQ.add_subspace_relation(OccupiedOrbital, InactiveOccupiedOrbital)

SASQ.add_space_intersection(ActiveOrbital, OccupiedOrbital, ActiveOccupiedOrbital)
SASQ.add_space_intersection(ActiveOrbital, VirtualOrbital, ActiveVirtualOrbital)
SASQ.add_space_intersection(InactiveOrbital, OccupiedOrbital, InactiveOccupiedOrbital)
SASQ.add_space_intersection(InactiveOrbital, VirtualOrbital, InactiveOccupiedOrbital)


