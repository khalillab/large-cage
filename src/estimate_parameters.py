#!/usr/bin/env python


import sys
import tempfile
import arviz as az
import numpy as np
import pymc3 as pm
import argparse
import subprocess
import pandas as pd


def get_options():
    description = ''
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('eggs',
                        help='File with empirical data')

    parser.add_argument('--cores',
                        type=int,
                        default=1,
                        help='Number of cores (default: %(default)d)')

    return parser.parse_args()


def load_observed(fname):
    edf = pd.read_csv(fname, sep='\t')

    edf = edf.iloc[1:].copy()
    edf = edf.rename(columns={'Unnamed: 1': 'day',
                              'Eggs': 'Cage A',
                              'Unnamed: 3': 'Cage B'})
    edf = edf[['day', 'Cage A', 'Cage B']].astype(int)
    data = np.concatenate([edf['Cage A'].values,
                           edf['Cage B'].values])
    return edf, data


if __name__ == '__main__':
    options = get_options()

    # here load actual eggs data
    edf, data = load_observed(options.eggs)

    # here define the simulator function
    # defined here so we can use the actual data
    # to parse the results
    def large_cage(a, b, edf=edf):
        tfile = tempfile.NamedTemporaryFile()

        command = '''python src/simulation.py \
        --drive 0 0 0 0 0 0 0 0 0 0 0 0 0 0 \
                0 0 0 0 0 0 0 0 0 0 0 0 0 0 \
                0 0 0 0 0 0 0 0 \
                228 228 228 228 228 228 \
        --antidote 0 0 0 0 0 0 0 0 0 0 0 0 0 0 \
                   0 0 0 0 0 0 0 0 0 0 0 0 0 0 \
                   0 0 0 0 0 0 0 0 \
                   0 0 0 0 0 0 \
        --wild-type 600 600 600 600 \
                    0 0 0 0 0 0 0 0 0 0 0 0 0 0 \
                    0 0 0 0 0 0 0 0 0 0 0 0 0 0 \
                    0 0 0 0 0 0 0 0 0 0 \
        --release 600 \
        --mating-probability %f \
        --egg-deposition-probability %f \
        --end-time 150 \
        --repetitions 2 > %s''' % (a, b, tfile.name)

        # run simulation
        subprocess.run(command, shell=True, check=True, stderr=subprocess.DEVNULL)
        
        # read output
        df = pd.read_csv(tfile.name, sep='\t')
        # time point 0 is GD release
        df['time'] = df['time'] - df[df['drives'] > 0]['time'].min()
        # select only the days for which we have data
        df = df[df['time'].isin(edf['day'].values)]
        # sort so that it matches the actual data
        data = df.sort_values(['round', 'time'])['eggs'].values

        return data
    
    with pm.Model() as simulated_cage:
        # here define parameters to estimate
        a = pm.Uniform('a', lower=0.1, upper=0.6)
        b = pm.Uniform('b', lower=0.05, upper=0.4)

        # here define how we run the simulation
        s = pm.Simulator('large-cage',
                         large_cage, params=(a, b),
                         epsilon=1,
                         observed=data)

        # actual sequential MC analysis
        trace, sim_data = pm.sample_smc(kernel='ABC',
                                        parallel=True,
                                        cores=options.cores,
                                        save_sim_data=True)
        idata = az.from_pymc3(trace, posterior_predictive=sim_data)

        summary = az.summary(idata, kind='stats')
        summary.to_csv(sys.stdout, sep='\t')
