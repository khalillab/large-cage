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
