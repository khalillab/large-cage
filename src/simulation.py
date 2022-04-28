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
    parser.add_argument('--special-release',
                        nargs='+',
                        default=None,
                        help='Release size for specific days; '
                             'format: day:release_size '
                             '(default: use --release-size)')
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
    parser.add_argument('--hom-antidote',
                        default=False,
                        action='store_true',
                        help='Release homozygous antidote males '
                             '(default: het.)')
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
    parser.add_argument('--mating-probability',
                        type=float,
                        default=agent.MATING_PROBABILITY,
                        help='Mating probability (default: %(default).2f)')
    parser.add_argument('--egg-deposition-probability',
                        type=float,
                        default=agent.EGG_DEPOSITION_PROBABILITY,
                        help='Egg deposition probability (default: %(default).2f)')
    parser.add_argument('--use-adults',
                        default=False,
                        action='store_true',
                        help='Use adults if no pupae are available '
                             '(useful for single releases, default: pupae only)')
    parser.add_argument('--end-time',
                        type=int,
                        default=365,
                        help='End day (default: %(default)d)')
    parser.add_argument('--time-points',
                        default=None,
                        help='Time points for printing simulation status '
                             '(one float per line; default: every day)')
    parser.add_argument('--late-releases',
                        type=int,
                        default=None,
                        help='Number of late releases (e.g. how many releases; '
                             'default: none, -1 indicates continuos releases)')
    parser.add_argument('--late-releases-start',
                        type=float,
                        default=30,
                        help='Start of late releases (e.g. what time; '
                             'default: %(default).2f, nothing happens unless '
                             'the --late-antidote or --late-wt args are used)')
    parser.add_argument('--late-antidote',
                        type=int,
                        default=None,
                        help='Late release of antidote adults '
                             '(how many adults to release; '
                             'default: no late release)')
    parser.add_argument('--late-wt',
                        type=int,
                        default=None,
                        help='Late release of wild-type adults '
                             '(how many adults to release; '
                             'default: no late release)')
    parser.add_argument('--eggs-filter',
                        default=False,
                        action='store_true',
                        help='Egg filter '
                             '(to reduce the number of eggs; '
                             'default: no filter)')
    parser.add_argument('--eggs-filter-mean',
                        type=float,
                        default=agent.EGGS_FILTER_MEAN,
                        help='Egg filter normal distribution mean '
                             '(to reduce the number of eggs; '
                             'default: %(default).3f)')
    parser.add_argument('--eggs-filter-std',
                        type=float,
                        default=agent.EGGS_FILTER_STD,
                        help='Egg filter normal distribution std-dev '
                             '(to reduce the number of eggs; '
                             'default: %(default).3f)')

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
    sys.stderr.write('Changing parameter EGGS_FILTER_MEAN from its '
                     'default (using the script arguments)\n')
    agent.EGGS_FILTER_MEAN = options.eggs_filter_mean
    sys.stderr.write('Changing parameter EGGS_FILTER_STD from its '
                     'default (using the script arguments)\n')
    agent.EGGS_FILTER_STD = options.eggs_filter_std
    sys.stderr.write('Changing parameter MATING_PROBABILITY from its '
                     'default (using the script arguments)\n')
    agent.MATING_PROBABILITY = options.mating_probability
    sys.stderr.write('Changing parameter EGG_DEPOSITION_PROBABILITY from its '
                     'default (using the script arguments)\n')
    agent.EGG_DEPOSITION_PROBABILITY = options.egg_deposition_probability

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

    special_releases = {}
    if options.special_release is not None:
        for release in options.special_release:
            day, size = release.split(':')
            special_releases[int(day)] = int(size)
    
    if options.eggs_filter:
        eggs_filter = (options.eggs_filter_mean,
                       options.eggs_filter_std)
    else:
        eggs_filter = None

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
                if not options.hom_antidote:
                    # het. male antidotes
                    x = Individual('m', ['W', 'W'], ['A', 'W'])
                else:
                    # hom. male antidotes
                    x = Individual('m', ['W', 'W'], ['A', 'A'])
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
       
        late_releases = None
        if options.late_releases is not None:
            late_releases_start = options.late_releases_start
            late_releases_counter = options.late_releases
            late = set()
            if options.late_antidote is not None:
                for i in range(int(options.late_antidote)):
                    if not options.hom_antidote:
                        # het. male antidotes
                        x = Individual('m', ['W', 'W'], ['A', 'W'])
                    else:
                        # hom. male antidotes
                        x = Individual('m', ['W', 'W'], ['A', 'A'])
                    x.stage = 'adult'
                    late.add(x)
            if options.late_wt is not None:
                for i in range(int(options.late_wt)):
                    x = Individual('m', ['W', 'W'], ['W', 'W'])
                    x.stage = 'adult'
                    late.add(x)
            late_releases = (late, late_releases_start, late_releases_counter)

        run_simulation(start_populations,
                       end_time=options.end_time,
                       repetition=j,
                       report_times=time_points,
                       release=agent.RELEASE,
                       special_releases=special_releases,
                       additional_releases=late_releases,
                       eggs_filter=eggs_filter,
                       time_step=agent.TIME_STEP,
                       release_days=agent.RELEASE_DAYS,
                       use_adults_if_needed=options.use_adults)
