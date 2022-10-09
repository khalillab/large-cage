#!/usr/bin/env python


import sys
import yaml
import argparse
from copy import deepcopy

from large_cage import agent
from large_cage.agent import get_all_genotypes
from large_cage.agent import Individual
from large_cage.agent import run_simulation
from large_cage.agent import print_header


def get_options():
    description = 'Run an individual-based simulation'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('--parameters',
                        default=None,
                        help='Override parameters using the ones found in '
                             'the provided YAML file. '
                             'NOTE: if a default parameter is not '
                             'in the provided file, then the default value '
                             'is used '
                             '(default: use the default ones)')

    return parser.parse_args()


if __name__ == "__main__":
    options = get_options()

    p = deepcopy(agent.params)

    # should we override the parameters?
    if options.parameters is not None:
        new_params = yaml.load(open(options.parameters),
                               Loader=yaml.SafeLoader)
        for k, v in new_params.items():
            if k in p:
                sys.stderr.write(f'Changing parameter {k} from its default\n')
            p[k] = v

    # check drive and antidote have the same length
    if len(p['RELEASE_WT']) != len(p['RELEASE_DRIVE']) or len(p['RELEASE_DRIVE']) != len(p['RELEASE_ANTI']):
        sys.stderr.write('Please provide the same number of introductions '
                         'for wild-type, drive and antidote individuals\n')
        sys.exit(1)
    WILDS = [int(x) for x in p['RELEASE_WT']]
    DRIVES = [int(x) for x in p['RELEASE_DRIVE']]
    ANTIDOTES = [int(x) for x in p['RELEASE_ANTI']]

    eggs_filter = None

    print_header(p)
    for j in range(p['REPETITIONS']):
        # init
        start_populations = []
        gd_populations = []
        for wild, drive, anti in zip(WILDS, DRIVES, ANTIDOTES):
            population = set()

            for i in range(drive):
                # het. male drives
                x = Individual('m', ['W', 'D'], ['W', 'W'],
                               nucl_from_father=True,
                               parameters=p)
                x.stage = 'adult'
                population.add(x)
            for i in range(anti):
                if not p['HOM_ANTIDOTE']:
                    # het. male antidotes
                    x = Individual('m', ['W', 'W'], ['A', 'W'],
                                   parameters=p)
                else:
                    # hom. male antidotes
                    x = Individual('m', ['W', 'W'], ['A', 'A'],
                                   parameters=p)
                x.stage = 'adult'
                population.add(x)

            # wild-type individuals
            for i in range(int(wild / 2)):
                x = Individual('f', ['W', 'W'], ['W', 'W'],
                               parameters=p)
                x.stage = 'adult'
                population.add(x)
                x = Individual('m', ['W', 'W'], ['W', 'W'],
                               parameters=p)
                x.stage = 'adult'
                population.add(x)

            start_populations.append(population)

        late_releases = None
        if p['LATE_RELEASES'] is not None:
            late_releases_start = p['LATE_RELEASES_START']
            late_releases_counter = p['LATE_RELEASES']
            late = set()
            if p['LATE_RELEASES_ANTI'] is not None:
                for i in range(int(p['LATE_RELEASES_ANTI'])):
                    if not p['HOM_ANTIDOTE']:
                        # het. male antidotes
                        x = Individual('m', ['W', 'W'], ['A', 'W'],
                                       parameters=p)
                    else:
                        # hom. male antidotes
                        x = Individual('m', ['W', 'W'], ['A', 'A'],
                                       parameters=p)
                    x.stage = 'adult'
                    late.add(x)
            late_releases = (late, late_releases_start, late_releases_counter)

        run_simulation(start_populations,
                       end_time=p['END_TIME'],
                       repetition=j,
                       report_times=None,
                       release=p['RELEASE'],
                       special_releases={},
                       additional_releases=late_releases,
                       eggs_filter=eggs_filter,
                       time_step=p['TIME_STEP'],
                       release_days=p['RELEASE_DAYS'],
                       use_adults_if_needed=p['USE_ADULTS'],
                       p=p)
