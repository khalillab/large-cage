#!/usr/bin/env python


import os
import sys
import numpy as np
import argparse
import pandas as pd
from sklearn import metrics


def get_options():
    description = ''
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('eggs',
                        help='File with empirical data')
    parser.add_argument('output',
                        help='Output file for an individual simulation')

    return parser.parse_args()


def load_observed(fname):
    edf = pd.read_csv(fname, sep='\t')

    edf = edf.iloc[1:].copy()
    edf = edf.rename(columns={'Unnamed: 1': 'day',
                              'Eggs': 'Cage A',
                              'Unnamed: 3': 'Cage B'})
    edf = edf[['day', 'Cage A', 'Cage B']].astype(int)
    # TODO: fix this
    # clip last time point, simulations did not include it
    edf = edf[edf['day'] <= 0].copy()
    data = np.concatenate([edf['Cage A'].values,
                           edf['Cage B'].values])
    return edf, data


def evaluate(values, data):
    if values.shape[0] * 2 < data.shape[0]:
        # simulation ended early
        # no score can be computed
        return np.nan

    y_hat = np.concatenate([values['eggs'].values,
                            values['eggs'].values])
    r2 = metrics.r2_score(data, y_hat)

    return r2


if __name__ == '__main__':
    options = get_options()

    # here load actual eggs data
    edf, real_data = load_observed(options.eggs)

    # read output
    df = pd.read_csv(options.output, sep='\t')
    # time point 0 is GD release
    df['time'] = df['time'] - df[df['drives'] > 0]['time'].min()
    # save average fitness value
    fitness = df[df['time'] < 0].groupby('round')['fitness'].mean().mean()
    # select only the days for which we have data
    df = df[df['time'].isin(edf['day'].values)]

    if df.shape[0] == 0:
        r2 = np.nan
    else:
        # compute R^2 for each round
        r2 = df.groupby('round').apply(evaluate, data=real_data)
        r2 = r2.mean()

    # derive parameters values
    mating, deposition = os.path.split(options.output)[-1].split('-')
    deposition = deposition.split('.tsv')[0]

    mating = float(mating)
    deposition = float(deposition)
    theoretical_fitness = mating * deposition

    print(f'{mating}\t{deposition}\t{theoretical_fitness}\t{fitness}\t{r2}')
