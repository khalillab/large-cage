#!/usr/bin/env python


import os
import yaml
import argparse
import numpy as np


def get_options():
    description = 'Generate parameter sweep YAML files'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('parameters',
                        help='YAML file to change')
    parser.add_argument('--out',
                        default='.',
                        help='Output directory (default: %(default)s)')

    return parser.parse_args()


if __name__ == "__main__":
    options = get_options()

    for p in range(500, 5100, 500):
        antidote = p * 0.23
        for d in np.linspace(0.1, 0.8, 15):
            y = yaml.load(open(options.parameters),
                      Loader=yaml.SafeLoader)

            y['RELEASE'] = p
            y['LATE_RELEASES_DRIVE_FREQUENCY'] = d
            fname = '.'.join(os.path.split(options.parameters)[-1].split('.')[:-1])
            try:
                os.mkdir(options.out)
            except: pass
            try:
                os.mkdir(os.path.join(options.out, 'release_%d_%.2f' % (p, d)))
            except: pass
            out = os.path.join(options.out, 'release_%d_%.2f' % (p, d), '%s.yaml' % fname)
            print(p, out)
            yaml.dump(y, open(out, 'w'))
