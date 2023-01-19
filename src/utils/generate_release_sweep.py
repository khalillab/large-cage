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
    parser.add_argument('--ignore-drive',
                        action='store_true',
                        default=False,
                        help='Don\'t sweep on drive frequency thresholds')
    parser.add_argument('--out',
                        default='.',
                        help='Output directory (default: %(default)s)')

    return parser.parse_args()


if __name__ == "__main__":
    options = get_options()

    drive_sweep = np.linspace(0.1, 0.8, 15)
    if options.ignore_drive:
        drive_sweep = [np.nan]

    for p in range(500, 5100, 500):
        antidote = int(p * 0.23)
        drive = int(p * 0.38)
        for d in drive_sweep:
            y = yaml.load(open(options.parameters),
                          Loader=yaml.SafeLoader)

            y['RELEASE'] = p
            y['RELEASE_WT'] = [p if x != 0
                               else 0
                               for x in y['RELEASE_WT']]
            y['RELEASE_DRIVE'] = [drive if x != 0
                                  else 0
                                  for x in y['RELEASE_DRIVE']]
            y['RELEASE_ANTI'] = [antidote if x != 0
                                 else 0
                                 for x in y['RELEASE_ANTI']]
            y['LATE_RELEASES_ANTI'] = antidote
            if not np.isnan(d):
                y['LATE_RELEASES_DRIVE_FREQUENCY'] = float(d)
            fname = '.'.join(os.path.split(options.parameters)[-1].split('.')[:-1])
            try:
                os.mkdir(options.out)
            except: pass
            try:
                if not np.isnan(d):
                    os.mkdir(os.path.join(options.out, 'release_%d_%.2f' % (p, d)))
                else:
                    os.mkdir(os.path.join(options.out, 'release_%d' % p))
            except: pass
            if not np.isnan(d):
                out = os.path.join(options.out, 'release_%d_%.2f' % (p, d), '%s.yaml' % fname)
            else:
                out = os.path.join(options.out, 'release_%d' % p, '%s.yaml' % fname)
            print(p, out)
            yaml.dump(y, open(out, 'w'))
