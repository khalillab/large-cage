#!/usr/bin/env python

import argparse

from large_cage.agent import REPETITIONS
from large_cage.agent import get_all_genotypes
from large_cage.agent import Individual
from large_cage.agent import run_simulation
from large_cage.agent import print_header 


def get_options():
    description = 'Run an individual-based simulation'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('--repetitions',
                        type=int,
                        default=REPETITIONS,
                        help='Repeated runs (default: %(default)d)')
    parser.add_argument('--pop-size',
                        type=int,
                        default=400,
                        help='Population size for initial release '
                             '(default: %(default)d)')
    parser.add_argument('--drive-1',
                        type=int,
                        default=71,
                        help='Population size for het. male drives, initial release '
                             '(default: %(default)d)')
    parser.add_argument('--drive-2',
                        type=int,
                        default=72,
                        help='Population size for het. male drives, second release '
                             '(default: %(default)d)')
    parser.add_argument('--time-points',
                        default=None,
                        help='Time points for printing simulation status '
                             '(one float per line; default: every day)')

    return parser.parse_args()

if __name__ == "__main__":
    options = get_options()

    # update parameters
    REPETITIONS = options.repetitions
    POPULATION = options.pop_size
    DRIVE_1 = options.drive_1
    DRIVE_2 = options.drive_2
    if options.time_points is not None:
        time_points = {float(l.rstrip()) for l in open(options.time_points)}
    else:
        time_points = None
    
    print_header()
    for j in range(REPETITIONS):
        # init
        start_populations = []
        gd_populations = []
        for drive in (DRIVE_1, DRIVE_2):
            population = set()

            for i in range(drive):
                x = Individual('m', ['W', 'D'], ['W', 'W'], 0, 0,
                               nucl_from_father=True)
                x.stage = 'adult'
                population.add(x)

            for i in range(int(POPULATION / 2)):
                x = Individual('f', ['W', 'W'], ['W', 'W'], 0, 0)
                x.stage = 'adult'
                population.add(x)
                x = Individual('m', ['W', 'W'], ['W', 'W'], 0, 0)
                x.stage = 'adult'
                population.add(x)

            start_populations.append(population)
        
        run_simulation(start_populations,
                       repetition=j,
                       report_times=time_points,
                       release=int(POPULATION / 2))
