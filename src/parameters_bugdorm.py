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
# days of the week for releases
# use integers here: 0 is Monday, six is Sunday
RELEASE_DAYS = [0,]
