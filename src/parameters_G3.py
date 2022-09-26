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
'f': {
     'DWAA': 0.82,
     'DWAW': 0.82,
     'DWWW': 0.45,
     'RWAW': 1.03,
     'WWAW': 1.03
     },
'm': {'DDAA': 0.1,
     'DRAA': 0.1,
     'DRAW': 1.06,
     'DRWW': 0.86,
     'DWAA': 0.1,
     'DWAW': 1.37,
     'DWWW': 1.28,
     'RRAA': 0.1,
     'RRAW': 0.94,
     'RWAA': 0.1,
     'RWAW': 0.94,
     'WWAA': 0.1,
     'WWAW': 0.81}
}

