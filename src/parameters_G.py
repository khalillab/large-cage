#!/usr/bin/env python

from scipy import stats as _stats
from large_cage.distributions import GaussianMixture as _GM

# scar code in case I misunderstood the scenario
#eggs_wt = [[_stats.norm(loc=90, scale=10),
#            _stats.norm(loc=1500, scale=300)],
#           [100000, int(100000 * 0.05)]]
#eggs_nucl_from_mother = [[_stats.norm(loc=78, scale=10),
#                          _stats.norm(loc=1299, scale=300)],
#                         [100000, int(100000 * 0.05)]]
#eggs_nucl_from_father = [[_stats.norm(loc=39, scale=10),
#                          _stats.norm(loc=651, scale=300)],
#                         [100000, int(100000 * 0.05)]]
#eggs_nucl_from_both = [[_stats.norm(loc=39, scale=10),
#                        _stats.norm(loc=651, scale=300)],
#                       [100000, int(100000 * 0.05)]]
#

_eggs_wt = [[_stats.norm(loc=0, scale=0),
             _stats.norm(loc=1500, scale=300)],
            [100000, int(100000 * 0.05)]]
_eggs_nucl_from_mother = [[_stats.norm(loc=0, scale=0),
                           _stats.norm(loc=1299, scale=300)],
                          [100000, int(100000 * 0.05)]]
_eggs_nucl_from_father = [[_stats.norm(loc=0, scale=0),
                           _stats.norm(loc=651, scale=300)],
                          [100000, int(100000 * 0.05)]]
_eggs_nucl_from_both = [[_stats.norm(loc=0, scale=0),
                         _stats.norm(loc=651, scale=300)],
                        [100000, int(100000 * 0.05)]]

# number of eggs
EGGS_WT = _GM(_eggs_wt[0],
              _eggs_wt[1])
EGGS_NUCL_FROM_MOTHER = _GM(_eggs_nucl_from_mother[0],
                            _eggs_nucl_from_mother[1])
EGGS_NUCL_FROM_FATHER = _GM(_eggs_nucl_from_father[0],
                            _eggs_nucl_from_father[1])
EGGS_NUCL_FROM_BOTH = _GM(_eggs_nucl_from_both[0],
                          _eggs_nucl_from_both[1])
