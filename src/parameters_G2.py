#!/usr/bin/env python

from scipy import stats as _stats

# parameters 
#
# NOTE: change values here and call the simulation.py script
#       with the --override-parameters src/parameters.py argument
#
# NOTE: if a parameter is removed from this file, the default in
#       src/large_cage/agent.py will be used
#
# genotype modifiers for mating probability
MATING_MOD = {
'f': {'DWAA': 0.82,
     'DWAW': 0.82,
     'DWWW': 0.45},
'm': {'DDAA': 0.69,
     'DDWW': 0.86,
     'DRAA': 0.69,
     'DRWW': 0.86,
     'DWAA': 0.69,
     'RRAA': 0.69,
     'RRAW': 0.94,
     'RWAA': 0.69,
     'RWAW': 0.94,
     'WWAA': 0.62,
     'WWAW': 0.81}
}
# genotype modifiers for number of eggs
EGGS_MOD = {
'f': {'DWAA': 0.96,
     'DWAW': 0.96,
     'DWWW': 0.65,},
'm': {}
}
# genotype modifiers for hatching rate
HATCHING_MOD = {
'f': {'DWAA': 0.81,
     'DWAW': 0.81,
     'DWWW': 0.62,},
'm': {}
}

