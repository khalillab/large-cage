#!/usr/bin/env python


import os
import yaml
import argparse


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

    y = yaml.load(open(options.parameters),
                  Loader=yaml.SafeLoader)

    releases = [600, 1200, 2400]
    times = [[1, 4], [1, 4, 7],
             [0, 1, 2, 3, 4, 5, 6]]

    for r in releases:
        for t in times:
            # this is the default scenario
            if r == 600 and len(t) == 2:
                continue
            y['RELEASE'] = r
            y['RELEASE_DAYS'] = t

            fname = '.'.join(os.path.split(options.parameters)[-1].split('.')[:-1])
            try:
                os.mkdir(options.out)
            except: pass
            try:
                os.mkdir(os.path.join(options.out, 'release_%d_%d' % (r, len(t))))
            except: pass
            out = os.path.join(options.out, 'release_%d_%d' % (r, len(t)), '%s.yaml' % fname)
            print(r, t, out)
            yaml.dump(y, open(out, 'w'))
