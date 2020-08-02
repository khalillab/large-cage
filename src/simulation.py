#!/usr/bin/env python


import sys
import argparse

from large_cage import agent
from large_cage.agent import get_all_genotypes
from large_cage.agent import Individual
from large_cage.agent import run_simulation
from large_cage.agent import print_header 


def get_options():
    description = 'Run an individual-based simulation'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('--override-parameters',
                        default=False,
                        action='store_true',
                        help='Override parameters using the ones found in '
                             'the parameters.py file that is in the same '
                             'directory as this script (default: use default ones)')
    parser.add_argument('--repetitions',
                        type=int,
                        default=agent.REPETITIONS,
                        help='Repeated runs (default: %(default)d)')
    parser.add_argument('--release-size',
                        type=int,
                        default=agent.RELEASE,
                        help='Release size (how many pupae to release; '
                             'default: %(default)d)')
    parser.add_argument('--wild-type',
                        nargs='+',
                        default=['400', '400'],
                        help='Population sizes for wild-type individuals, initial releases '
                             '(each initial introduction; half are males, half females, '
                             'default: %(default)s)')
    parser.add_argument('--drive',
                        nargs='+',
                        default=['71', '72'],
                        help='Population sizes for het. male drives, initial releases '
                             '(each initial introduction; default: %(default)s)')
    parser.add_argument('--antidote',
                        nargs='+',
                        default=['0', '0'],
                        help='Population sizes for het. male antidotes, initial releases '
                             '(each initial introduction; default: %(default)s)')
    parser.add_argument('--hom-antidote-effect',
                        type=float,
                        default=agent.HOM_ANTIDRIVE_EFFECT,
                        help='Eggs output reduction for hom. antidotes '
                             '(default: %(default).2f)')
    parser.add_argument('--het-antidote-effect',
                        type=float,
                        default=agent.HET_ANTIDRIVE_EFFECT,
                        help='Eggs output reduction for het. antidotes '
                             '(default: %(default).2f)')
    parser.add_argument('--time-points',
                        default=None,
                        help='Time points for printing simulation status '
                             '(one float per line; default: every day)')

    return parser.parse_args()


if __name__ == "__main__":
    options = get_options()

    # update parameters
    sys.stderr.write('Changing parameter REPETITION from its default (using '
                     'the script arguments)\n')
    agent.REPETITIONS = options.repetitions
    sys.stderr.write('Changing parameter RELEASE from its default (using '
                     'the script arguments)\n')
    agent.RELEASE = options.release_size
    sys.stderr.write('Changing parameter HOM_ANTIDRIVE_EFFECT from its '
                     'default (using the script arguments)\n')
    agent.HOM_ANTIDRIVE_EFFECT = options.hom_antidote_effect
    sys.stderr.write('Changing parameter HET_ANTIDRIVE_EFFECT from its '
                     'default (using the script arguments)\n')
    agent.HET_ANTIDRIVE_EFFECT = options.het_antidote_effect

    # should we override the parameters?
    if options.override_parameters:
        try:
            import parameters
        except ImportError:
            sys.stderr.write('The parameters.py file should be in the same '
                             'directory as this script in order to override '
                             'default parameters\n')
            sys.exit(1)
        for var in dir(parameters):
            if var.startswith('_'):
                continue
            sys.stderr.write(f'Changing parameter {var} from its default\n')
            setattr(agent, var, getattr(parameters, var))

    # check drive and antidote have the same length
    if len(options.wild_type) != len(options.drive) or len(options.drive) != len(options.antidote):
        sys.stderr.write('Please provide the same number of introductions '
                         'for wild-type, drive and antidote individuals\n')
        sys.exit(1)
    WILDS = [int(x) for x in options.wild_type]
    DRIVES = [int(x) for x in options.drive]
    ANTIDOTES = [int(x) for x in options.antidote]
    if options.time_points is not None:
        time_points = {float(l.rstrip()) for l in open(options.time_points)}
    else:
        time_points = None
    
    print_header()
    for j in range(agent.REPETITIONS):
        # init
        start_populations = []
        gd_populations = []
        for wild, drive, anti in zip(WILDS, DRIVES, ANTIDOTES):
            population = set()

            for i in range(drive):
                # het. male drives
                x = Individual('m', ['W', 'D'], ['W', 'W'], 
                               nucl_from_father=True)
                x.stage = 'adult'
                population.add(x)
            for i in range(anti):
                # het. male antidotes
                x = Individual('m', ['W', 'W'], ['A', 'W'])
                x.stage = 'adult'
                population.add(x)

            # wild-type individuals
            for i in range(int(wild / 2)):
                x = Individual('f', ['W', 'W'], ['W', 'W'],)
                x.stage = 'adult'
                population.add(x)
                x = Individual('m', ['W', 'W'], ['W', 'W'],)
                x.stage = 'adult'
                population.add(x)

            start_populations.append(population)
        
        run_simulation(start_populations,
                       repetition=j,
                       report_times=time_points,
                       release=agent.RELEASE)
