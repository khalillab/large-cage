#!/usr/bin/env python

import numpy as np
from scipy import stats

class GaussianMixture():
    def __init__(self, distributions, samples):
        self._samples = np.concatenate([dist.rvs(n)
                                        for dist, n in zip(distributions,
                                                           samples)])
        self.kde = stats.gaussian_kde(self._samples)

    def rvs(self, n=1):
        samples = self.kde.resample(n)[0]
        if n == 1:
            return np.abs(samples[0])
        else:
            return np.abs(samples)
