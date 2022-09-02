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
'f': {'DDAA': 0.1,
     'DDAW': 0.1,
     'DDWW': 0.1,
     'DRAA': 0.1,
     'DWAA': 0.1,
     'DWAW': 0.64,
     'DWWW': 0.13,
     'RRAA': 0.1,
     'RWAA': 0.1,
     'RWAW': 1.19,
     'WWAA': 0.1,
     'WWAW': 1.19},
'm': {'DDAA': 0.1,
     'DDAW': 0.1,
     'DDWW': 0.1,
     'DRAA': 0.1,
     'DRAW': 1.06,
     'DRWW': 0.86,
     'DWAA': 0.1,
     'DWAW': 1.06,
     'DWWW': 0.86,
     'RRAA': 0.1,
     'RRAW': 0.94,
     'RWAA': 0.1,
     'RWAW': 0.94,
     'WWAA': 0.1,
     'WWAW': 0.94}
}

