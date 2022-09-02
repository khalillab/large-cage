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
'f': {'DWAA': 0.64,
     'DWAW': 0.64,
     'DWWW': 0.13},
'm': {'DDAA': 0.69,
     'DDWW': 0.86,
     'DRAA': 0.69,
     'DRWW': 0.86,
     'DWAA': 0.69,
     'DWWW': 0.86,
     'RRAA': 0.69,
     'RRAW': 0.94,
     'RWAA': 0.69,
     'RWAW': 0.94,
     'WWAA': 0.69,
     'WWAW': 0.94}
}
# genotype modifiers for number of eggs
EGGS_MOD = {
'f': {'DWAA': 0.96,
     'DWAW': 0.96,
     'DWWW': 0.65,},
'm': {}
}
# genotype modifiers for deposition probability
DEPOSITION_MOD = {
'f': {'DWAA': 0.58,
     'DWAW': 0.58,
     'DWWW': 0.5,
     'RWAA': 0.58,
     'RWAW': 0.58,
     'RWWW': 0.5,
     'WWAA': 0.9,
     'WWWW': 0.95},
'm': {'DDAA': 0.79,
     'DDAW': 0.89,
     'DDWW': 0.94,
     'DRAA': 0.79,
     'DRAW': 0.89,
     'DRWW': 0.94,
     'DWAA': 0.79,
     'DWAW': 0.89,
     'DWWW': 0.94,
     'RRAA': 0.79,
     'RRAW': 0.89,
     'RRWW': 0.94,
     'RWAA': 0.79,
     'RWAW': 0.89,
     'RWWW': 0.94,
     'WWAA': 0.79,
     'WWAW': 0.97,
     'WWWW': 0.95}
}
# genotype modifiers for hatching rate
HATCHING_MOD = {
'f': {'DWAA': 0.81,
     'DWAW': 0.81,
     'DWWW': 0.62,},
'm': {}
}
# lifespan for adults (weibull)
# a random variable is derived from this distribution
SURVIVAL = _stats.weibull_min(c=4.238761472584551,
                              scale=27.901772767965426,
                              loc=-5.4644363721937825)
SURVIVAL_MALE = _stats.weibull_min(c=1.709696e+00,
                                   scale=1.142531e+01,
                                   loc=1.833875e+00)
SURVIVAL_FEMALE = _stats.weibull_min(c=8.215818e+00,
                                     scale=4.835784e+01,
                                     loc=-2.551209e+01)
# days of the week for releases
# use integers here: 0 is Monday, six is Sunday
RELEASE_DAYS = [0,]
