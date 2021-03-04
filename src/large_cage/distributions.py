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

class DistributionMixture():
    def __init__(self, distributions, weights):
        self._distributions = distributions
        self._weights = np.array(weights)
        self._weights /= self._weights.sum()

    def rvs(self, n=1):
        samples = []
        for _ in range(n):
            dist = np.random.choice(np.arange(len(self._distributions)),
                                    p=self._weights)
            sample = np.abs(self._distributions[dist].rvs(1))
            samples.append(sample)
        if n == 1:
            return samples[0]
        else:
            return np.array(samples)
