#!/usr/bin/env python


import os
import yaml
import argparse


def get_options():
    description = 'Generate parameter sweep YAML files'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('parameters',
                        help='YAML file to change')
    parser.add_argument('--het',
                        action='store_true',
                        default=False,
                        help='Apply cost to heterozygous as well')
    parser.add_argument('--het-model',
                        choices=['epistasis+',
                                 'additive', 'same'],
                        default='same',
                        help='How to apply mating cost to antidote hets (%(default)s)')
    parser.add_argument('--transhet',
                        action='store_true',
                        default=False,
                        help='Do not apply cost to transhet individuals')
    parser.add_argument('--modify',
                        action='store_true',
                        default=False,
                        help='Modify existing cost (default: replace it)')
    parser.add_argument('--out',
                        default='.',
                        help='Output directory (default: %(default)s)')

    return parser.parse_args()


if __name__ == "__main__":
    options = get_options()

    sexes = ['m']
    locus1 = ['D', 'W', 'R']
    genotypes1 = {tuple(sorted((x, y)))
                  for x in locus1
                  for y in locus1}
    locus2 = ['A', 'W']
    genotypes2 = {tuple(sorted((x, y)))
                  for x in locus2
                  for y in locus2}
    genotypes = {''.join((x, y, w, z))
                 for x, y in genotypes1
                 for w, z in genotypes2}

    transhet = {x for x in genotypes
                if 'D' in x[:2]
                and 'A' in x[2:]}

    y = yaml.load(open(options.parameters),
                  Loader=yaml.SafeLoader)

    for p in [0.01, 0.05, 0.2, 0.3, 0.5,
              0.7, 0.9, 1]:
        if options.het_model == 'same':
            het_value = p
        elif options.het_model == 'additive':
            het_value = p * 2
        elif options.het_model == 'epistasis+':
            het_value = p * 4
            if het_value > 1:
                het_value = 1

        for sex in y['MATING_MOD']:
            if sex not in sexes:
                continue
            for genotype in genotypes:
                if options.transhet and genotype in transhet:
                    continue
                v = y['MATING_MOD'][sex].get(genotype, 1.)
                if genotype.endswith('AA'):
                    if options.modify:
                        y['MATING_MOD'][sex][genotype] = p * v
                    else:
                        y['MATING_MOD'][sex][genotype] = p

                elif 'A' in genotype[2:] and options.het:
                    if options.modify:
                        y['MATING_MOD'][sex][genotype] = het_value * v
                    else:
                        y['MATING_MOD'][sex][genotype] = het_value
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
