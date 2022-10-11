#!/usr/bin/env python


import argparse
import pandas as pd


def get_options():
    description = 'Convert tsv output to an excel multisheet file'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('input',
                        help='Input file')
    parser.add_argument('output',
                        help='Output file')

    return parser.parse_args()


if __name__ == "__main__":
    options = get_options()

    m = pd.read_csv(options.input, sep='\t')
    ex = pd.ExcelWriter(options.output)

    cols = ['round', 'time', 'pop', 'fpop', 'eggs', 'feggs', 'output',
            'WT.output', 'transgenes.output', 'drives.output',
            'anti.output', 'resistance.output']

    t = m[m['initial_release']][cols].pivot_table(index='time', columns='round')

    for c in cols:
        if c in ['round', 'time']:
            continue
        t[c].to_excel(ex, sheet_name=c)

    ex.close()
