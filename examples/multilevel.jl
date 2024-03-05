using SpinAdaptedSecondQuantization

# Define spaces with names and colors

ActiveOrbital = SASQ.new_space(:ActiveOrbital, "a", "pqrstuv")
InactiveOrbital = SASQ.new_space(:InactiveOrbital, "i", "PQRSTUV")

# set_color(ActiveOrbital, :cyan)
# set_color(InactiveOrbital, :red)

ActiveVirtualOrbital = SASQ.new_space(:ActiveVirtualOrbital, "av", "abcdefg")
ActiveOccupiedOrbital = SASQ.new_space(:ActiveOccupiedOrbital, "ao", "ijklmno")

# set_color(ActiveVirtualOrbital, :cyan)
# set_color(ActiveOccupiedOrbital, :cyan)

InactiveVirtualOrbital = SASQ.new_space(:InactiveVirtualOrbital, "iv", "ABCDEFG")
InactiveOccupiedOrbital = SASQ.new_space(:InactiveOccupiedOrbital, "io", "IJKLMNO")

# set_color(InactiveVirtualOrbital, :red)
# set_color(InactiveOccupiedOrbital, :red)

# Define relations

# SASQ.add_space_sum(ActiveOrbital, InactiveOrbital, GeneralOrbital)
SASQ.add_subspace_relation(GeneralOrbital, ActiveOrbital)
# SASQ.add_subspace_relation(GeneralOrbital, InactiveOrbital)
SASQ.add_subspace_relation(OccupiedOrbital, ActiveOccupiedOrbital)
SASQ.add_subspace_relation(VirtualOrbital, InactiveOrbital)
# SASQ.add_subspace_relation(VirtualOrbital, InactiveVirtualOrbital)

# SASQ.add_space_sum(ActievVirtualOrbital, InactiveVirtualOrbital, VirtualOrbital)
# SASQ.add_space_sum(ActiveOccupiedOrbital, InactiveOccupiedOrbital, OccupiedOrbital)
SASQ.add_space_sum(ActiveOccupiedOrbital, ActiveVirtualOrbital, ActiveOrbital)
SASQ.add_space_sum(InactiveOccupiedOrbital, InactiveVirtualOrbital, InactiveOrbital)

SASQ.add_space_intersection(ActiveOrbital, OccupiedOrbital, ActiveOccupiedOrbital)
SASQ.add_space_intersection(ActiveOrbital, VirtualOrbital, ActiveVirtualOrbital)
# SASQ.add_space_intersection(InactiveOrbital, OccupiedOrbital, InactiveOccupiedOrbital)
# SASQ.add_space_intersection(InactiveOrbital, VirtualOrbital, InactiveVirtualOrbital)

# Macros for constrain

active(indices...) = constrain(p => ActiveOrbital for p in indices)
inactive(indices...) = constrain(p => InactiveOrbital for p in indices)

aocc(indices...) = constrain(p => ActiveOccupiedOrbital for p in indices)
iocc(indices...) = constrain(p => InactiveOccupiedOrbital for p in indices)
avir(indices...) = constrain(p => ActiveVirtualOrbital for p in indices)
ivir(indices...) = constrain(p => InactiveVirtualOrbital for p in indices)
