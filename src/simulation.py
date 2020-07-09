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

    parser.add_argument('--repetitions',
                        type=int,
                        default=agent.REPETITIONS,
                        help='Repeated runs (default: %(default)d)')
    parser.add_argument('--pop-size',
                        type=int,
                        default=400,
                        help='Population size for initial release '
                             '(default: %(default)d)')
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
    parser.add_argument('--time-points',
                        default=None,
                        help='Time points for printing simulation status '
                             '(one float per line; default: every day)')

    return parser.parse_args()


if __name__ == "__main__":
    options = get_options()

    # update parameters
    agent.REPETITIONS = options.repetitions
    agent.POPULATION = options.pop_size
    
    # check drive and antidote have the same length
    if len(options.drive) != len(options.antidote):
        sys.stderr.write('Please provide the same number of introductions '
                         'for drive and antidote individuals\n')
        sys.exit(1)
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
        for drive, anti in zip(DRIVES, ANTIDOTES):
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
            for i in range(int(agent.POPULATION / 2)):
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
                       release=agent.POPULATION)
