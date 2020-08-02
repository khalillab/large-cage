#!/usr/bin/env python


import sys
import argparse
import pandas as pd


def get_options():
    description = 'Give average and MAD for multiple simulations'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('table',
                        help='Model\'s output (tsv format)')

    return parser.parse_args()


if __name__ == "__main__":
    options = get_options()

    m = pd.read_csv(options.table, sep='\t')
    avg = m.groupby(['time']).mean().drop(columns=['round'])
    mad = m.groupby(['time']).mean().drop(columns=['round', 'initial_release'])
    mad.columns = [f'{x}.MAD' for x in mad.columns]
    m = avg.join(mad, how='outer').reset_index()
    m.to_csv(sys.stdout, sep='\t', index=False)
