#!/usr/bin/env python


import argparse
import pandas as pd


def get_options():
    description = 'Merge multiple runs in a single file'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('input',
                        nargs='+',
                        help='Input file')
    parser.add_argument('output',
                        help='Output file')

    return parser.parse_args()


if __name__ == "__main__":
    options = get_options()

    res = []
    for i, f in enumerate(options.input):
        m = pd.read_csv(f, sep='\t')
        m['round'] = i
        res.append(m)

    pd.concat(res).to_csv(options.output, sep='\t', index=False)
