#!/usr/bin/env python


import sys
import argparse
import pandas as pd


def get_options():
    description = 'Select only relevant columns from model\'s output'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('table',
                        help='Model\'s output (tsv format)')

    return parser.parse_args()


if __name__ == "__main__":
    options = get_options()

    m = pd.read_csv(options.table, sep='\t')
    m = m[['round', 'time', 'initial_release', 'pop', 'fpop', 'eggs',
           'feggs', 'output', 'foutput', 
           'WT.output', 'transgenes.output', 'drives.output', 'anti.output',
           'resistance.output', 'DDAA.output', 'DDAW.output', 'DDWW.output',
           'DRAA.output', 'DRAW.output', 'DRWW.output', 'DWAA.output',
           'DWAW.output', 'DWWW.output', 'RRAA.output', 'RRAW.output',
           'RRWW.output', 'RWAA.output', 'RWAW.output', 'RWWW.output',
           'WWAA.output', 'WWAW.output', 'WWWW.output']]
    m = m.rename(columns={x: x.replace('.output', '')
                          for x in m.columns
                          if x.endswith('.output')})
    m.to_csv(sys.stdout, sep='\t', index=False)
