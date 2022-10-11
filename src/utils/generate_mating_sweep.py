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

    for p in [0.01, 0.05, 0.2, 0.3, 0.5,
              0.7, 0.9, 1]:
        for sex in y['MATING_MOD']:
            for k, v in y['MATING_MOD'][sex].items():
                if k.endswith('AA'):
                    y['MATING_MOD'][sex][k] = p
        fname = '.'.join(os.path.split(options.parameters)[-1].split('.')[:-1])
        try:
            os.mkdir(options.out)
        except: pass
        try:
            os.mkdir(os.path.join(options.out, 'lowfitness_%.2f' % p))
        except: pass
        out = os.path.join(options.out, 'lowfitness_%.2f' % p, '%s.yaml' % fname)
        print(p, out)
        yaml.dump(y, open(out, 'w'))
