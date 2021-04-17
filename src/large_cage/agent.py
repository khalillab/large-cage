#!/usr/bin/env python

import sys
import random
import itertools
import numpy as np
from scipy import stats
from copy import deepcopy


# parameters 
#
# can be overridden by calling scipts
#
# release size
# (i.e. how many pupae)
RELEASE = 400
# how many simulations to run
REPETITIONS = 100
# fitness effect of the antidrive (on number of eggs)
HOM_ANTIDRIVE_EFFECT = 0.15
HET_ANTIDRIVE_EFFECT = 0.3
# females mating probability
MATING_PROBABILITY = 0.6
# can females/males mate multiple times?
MULTIPLE_MATING_FEMALE = False
MULTIPLE_MATING_MALE = True
# females eggs deposition probability
EGG_DEPOSITION_PROBABILITY = 0.5
# lifespan for adults (weibull)
# a random variable is derived from this distribution
SURVIVAL = stats.weibull_min(c=2.2472084592310644,
                             scale=6.213064151445494,
                             loc=0.7275226197070571)
# times from egg to full adult
# time 0 is egg deposition
# timer is reset once adult stage is reached
# times are extracted from a uniform distribution
# between the two extremes
TIME_TO_HATCH = (1.8, 2.2)
TIME_TO_PUPA = (9, 10)
TIME_TO_MATURATION = (10, 12)
# frequency of update of the simulation
# best not to change this in order to avoid missing
# introductions and other unforeseen bugs
TIME_STEP = 0.1
# days of the week for releases
# use integers here: 0 is Monday, six is Sunday
RELEASE_DAYS = [1, 4]
# alleles (at two loci)
WILD_TYPE = 'W' # both loci
DRIVE = 'D' # locus 1
ANTI_DRIVE = 'A' # locus 2
RESISTANCE = 'R' # locus 1
# which genotypes are sterile (locus 1)
NON_FUNCTIONAL = [
        {'D'}, # homozygous
        {'D', 'R'}, # heterozygous
        {'R'}] # homozygous
# phenotype of different genotypes
# for each parameter we define the normal distribution
# from which to extract the random variable
# 
# NUCL_FROM_* indicates from which parent the nuclease
# has been inherited
#
# homing efficiency
DRIVE_EFFICIENCY_NUCL_FROM_FATHER = stats.norm(loc=0.9867,
                                               scale=0.0163)
DRIVE_EFFICIENCY_NUCL_FROM_MOTHER = stats.norm(loc=0.9667,
                                               scale=0.0408)
DRIVE_EFFICIENCY_NUCL_FROM_BOTH = 1.
# antidote efficiency
ANTI_DRIVE_EFFICIENCY = {'f': 1.,
                         'm': 1.}
# emergence of resistant allele (if drive fails)
RESISTANCE_EFFICIENCY = {'f': 0.4685,
                         'm': 0.4685}
# number of eggs
EGGS_WT = stats.norm(loc=137.4,
                     scale=34.5)
EGGS_NUCL_FROM_MOTHER = stats.norm(loc=118.96,
                                   scale=34.5)
EGGS_NUCL_FROM_FATHER = stats.norm(loc=59.67,
                                   scale=50.6)
EGGS_NUCL_FROM_BOTH = stats.norm(loc=59.67,
                                 scale=50.6)
# eggs hatching probability
HATCHING_WT = stats.norm(loc=0.8667,
                         scale=0.0046)
HATCHING_NUCL_FROM_MOTHER = stats.norm(loc=0.5313,
                                       scale=0.0168)
HATCHING_NUCL_FROM_FATHER = stats.norm(loc=0.8725,
                                       scale=0.0159)
HATCHING_NUCL_FROM_BOTH = stats.norm(loc=0.5094,
                                     scale=0.0553)
# larval mortality 
LARVAL_WT = stats.norm(loc=0.0825,
                       scale=0.0214)
LARVAL_NUCL_FROM_MOTHER = stats.norm(loc=0.1019,
                                     scale=0.0168)
LARVAL_NUCL_FROM_FATHER = stats.norm(loc=0.0671,
                                     scale=0.0146)
LARVAL_NUCL_FROM_BOTH = stats.norm(loc=0.0949,
                                   scale=0.0214)
# pupal mortality (males)
PUPAL_M_WT = stats.norm(loc=0.1837,
                        scale=0.0627)
PUPAL_M_NUCL_FROM_MOTHER = stats.norm(loc=0.0783,
                                      scale=0.0167)
PUPAL_M_NUCL_FROM_FATHER = stats.norm(loc=0.0692,
                                      scale=0.0233)
PUPAL_M_NUCL_FROM_BOTH = stats.norm(loc=0.0796,
                                    scale=0.0182)
# pupal mortality (females)
PUPAL_F_WT = stats.norm(loc=0.0918,
                        scale=0.0161)
PUPAL_F_NUCL_FROM_MOTHER = stats.norm(loc=0.0503,
                                      scale=0.0121)
PUPAL_F_NUCL_FROM_FATHER = stats.norm(loc=0.0676,
                                      scale=0.0223)
# is this individual intersex?
# (means that it does not mate)
# requires presence of the drive in locus 1
INTERSEX_NUCL_FROM_FATHER = stats.norm(loc=0.0353,
                                       scale=0.0124)
INTERSEX_NUCL_FROM_MOTHER = stats.norm(loc=0.0096,
                                       scale=0.0066)
# eggs filter
# (to simulate any problems with eggs production)
EGGS_FILTER_MEAN = 1750 
EGGS_FILTER_STD = 150

class Individual():
    '''Individual mosquito model

    Two loci are modelled in this simulation:
    - locus 1: dsx
    - locus 2: presence/absence of antidote

    Individuals are in the egg stage upon object creation,
    and the time to reach each stage is predetermined by drawing
    random variables from predefined distributions
    
    Mating, egg deposition probabilities, number of eggs produced
    and similar quantities are also drawn from predefined distributions

    Example 1: wild-type male
    >>> Individual('m', ['W', 'W'], ['W', 'W'])
    
    Example 2: het. drive male (nuclease received from father)
    >>> Individual('m', ['D', 'W'], ['W', 'W'], nucl_from_father=True)
    
    Example 3: hom. drive female (nuclease received from both parents)
    >>> Individual('f', ['D', 'D'], ['W', 'W'],
                   nucl_from_father=True,
                   nucl_from_mother=True)
    
    Example 4: het. drive/resistant female, het. antidrive (nuclease received from mother)
    >>> Individual('f', ['D', 'R'], ['A', 'W'],
                   nucl_from_mother=True)

    Useful methods:

    Increase the age of the individual by 1 day. Change of stage is automatically
    set by the method (including death)
    >>> i.change_age(1)

    Is this individual still alive?
    >>> i.is_alive()
    False

    Is this individual going to develop into an adult?
    >>> i.will_develop()
    True

    Form gamete for locus 1 (dsx)
    >>> i.form_gamete1()
    D

    Form gamete for locus 2 (antidote)
    >>> i.form_gamete2()
    W
    '''
    def __init__(self, sex, genotype1, genotype2,
                 nucl_from_father=False, nucl_from_mother=False):
        '''Create a new individual for the simulation
        Starts from the egg stage

        Args:
            sex (str)
                Sex: one of ['f', 'm']
            genotype1 (iterable)
                All unique alleles at the dsx locus
            genotype2 (iterable)
                All unique alleles at the antidote locus
            nucl_from_father (bool)
                Wether the father passes the nuclease to the egg
            nucl_from_mother (bool)
                Wether the mother passes the nuclease to the egg
        '''
        self.sex = sex
        self.nucl_from_father = nucl_from_father
        self.nucl_from_mother = nucl_from_mother

        if len(set(genotype1)) == 1:
            self.hom1 = True
        else:
            self.hom1 = False
        self.genotype1 = set(genotype1)
        if len(set(genotype2)) == 1:
            self.hom2 = True
        else:
            self.hom2 = False
        self.genotype2 = set(genotype2)

        self.age = 0.
        self.time_to_hatch = random.uniform(TIME_TO_HATCH[0], TIME_TO_HATCH[1])
        self.time_to_pupa = random.uniform(TIME_TO_PUPA[0], TIME_TO_PUPA[1])
        self.time_to_maturation = random.uniform(TIME_TO_MATURATION[0], TIME_TO_MATURATION[1])
        # lifespan after full maturation
        self.death = SURVIVAL.rvs()
        # egg -> larva -> pupa -> adult
        self.stage = 'egg'

        # keep track of mating events
        self.mated = False

        # is this individual intersex
        self.intersex = self._is_intersex()
        
        # initial mating probability 
        self.mating = self.get_mating()

        # will deposit eggs?
        self.deposing_eggs = self._eggs_deposition()
        
        # egg production
        if self.sex == 'f' and self.genotype1 not in NON_FUNCTIONAL:
            if DRIVE not in self.genotype1:
                self.eggs = int(EGGS_WT.rvs())
            elif DRIVE in self.genotype1 and nucl_from_father and nucl_from_mother:
                self.eggs = int(EGGS_NUCL_FROM_BOTH.rvs())
            elif DRIVE in self.genotype1 and nucl_from_father:
                self.eggs = int(EGGS_NUCL_FROM_FATHER.rvs())
            elif DRIVE in self.genotype1 and nucl_from_mother:
                self.eggs = int(EGGS_NUCL_FROM_MOTHER.rvs())
            else:
                raise RuntimeError(self.genotype1, nucl_from_father, nucl_from_mother)
        else:
            self.eggs = 0
        if ANTI_DRIVE in self.genotype2:
            if self.hom2:
                self.eggs = int(round(self.eggs * HOM_ANTIDRIVE_EFFECT))
            else:
                self.eggs = int(round(self.eggs * HET_ANTIDRIVE_EFFECT))

        # hatching probability
        if DRIVE not in self.genotype1:
            hatching = HATCHING_WT.rvs()
        elif DRIVE in self.genotype1 and nucl_from_father and nucl_from_mother:
            hatching = HATCHING_NUCL_FROM_BOTH.rvs()
        elif DRIVE in self.genotype1 and nucl_from_father:
            hatching = HATCHING_NUCL_FROM_FATHER.rvs()
        elif DRIVE in self.genotype1 and nucl_from_mother:
            hatching = HATCHING_NUCL_FROM_MOTHER.rvs()
        else:
            raise RuntimeError(self.genotype1, nucl_from_father, nucl_from_mother)
        if random.random() < hatching:
            self.hatching = True
        else:
            self.hatching = False

        # larval mortality
        if DRIVE not in self.genotype1:
            larval = LARVAL_WT.rvs()
        elif DRIVE in self.genotype1 and nucl_from_father and nucl_from_mother:
            larval = LARVAL_NUCL_FROM_BOTH.rvs()
        elif DRIVE in self.genotype1 and nucl_from_father:
            larval = LARVAL_NUCL_FROM_FATHER.rvs()
        elif DRIVE in self.genotype1 and nucl_from_mother:
            larval = LARVAL_NUCL_FROM_MOTHER.rvs()
        else:
            raise RuntimeError(self.genotype1, nucl_from_father, nucl_from_mother)
        if random.random() < larval:
            self.larva = False
        else:
            self.larva = True

        # pupal mortality
        if sex == 'm':
            if DRIVE not in self.genotype1:
                pupal = PUPAL_M_WT.rvs()
            elif DRIVE in self.genotype1 and nucl_from_father and nucl_from_mother:
                pupal = PUPAL_M_NUCL_FROM_BOTH.rvs()
            elif DRIVE in self.genotype1 and nucl_from_father:
                pupal = PUPAL_M_NUCL_FROM_FATHER.rvs()
            elif DRIVE in self.genotype1 and nucl_from_mother:
                pupal = PUPAL_M_NUCL_FROM_MOTHER.rvs()
            else:
                raise RuntimeError(self.genotype1, nucl_from_father, nucl_from_mother)
        else:
            if DRIVE not in self.genotype1:
                pupal = PUPAL_F_WT.rvs()
            elif DRIVE in self.genotype1 and nucl_from_father and nucl_from_mother:
                # assumed
                pupal = PUPAL_M_NUCL_FROM_BOTH.rvs()
            elif DRIVE in self.genotype1 and nucl_from_father:
                pupal = PUPAL_F_NUCL_FROM_FATHER.rvs()
            elif DRIVE in self.genotype1 and nucl_from_mother:
                pupal = PUPAL_F_NUCL_FROM_MOTHER.rvs()
            else:
                raise RuntimeError(self.genotype1, nucl_from_father, nucl_from_mother)
        if random.random() < pupal:
            self.pupa = False
        else:
            self.pupa = True

    def get_genotype(self):
        '''Get the full genotype at both loci

        Returns
           genotype (str)
               The first two chars are the dsx locus, the last the antidote

        Examples
        >>> i.get_genotype() # wild-type
        WWWW
        >>> i.get_genotype() # hom. drive, het. antidote
        DDAW
        >>> i.get_genotype() # het. drive, hom. antidote
        DWAA
        '''
        gt = ''
        for genotype in (self.genotype1, self.genotype2):
            if len(genotype) == 1:
                gt += list(genotype)[0]
                gt += list(genotype)[0]
            else:
                gt += ''.join(sorted(genotype))
        return gt

    def _is_intersex(self):
        if DRIVE in self.genotype1 and (self.nucl_from_father and self.nucl_from_mother):
            intersex = 0
        elif DRIVE in self.genotype1 and self.nucl_from_father and self.sex == 'f':
            intersex = INTERSEX_NUCL_FROM_FATHER.rvs()
        elif DRIVE in self.genotype1 and self.nucl_from_mother and self.sex == 'f':
            intersex = INTERSEX_NUCL_FROM_MOTHER.rvs()
        else:
            intersex = 0
        if random.random() <= intersex:
            return True
        return False

    def _mating_probability(self):
        return MATING_PROBABILITY

    def get_mating(self):
        '''Generate the probability that this individual will mate
        
        Returns
            mating (bool)
               Wether the individual will mate

        The decision is based on sex, genotype and intersex phenotype
        '''
        if self.intersex:
            return False
        if self.sex == 'm':
            return True
        if self.genotype1 in NON_FUNCTIONAL:
            return False
        else:
            prob = self._mating_probability()
            if random.random() < prob:
                return True
            else:
                return False

    def _eggs_deposition(self):
        if self.sex == 'm':
            return None
        elif random.random() <= EGG_DEPOSITION_PROBABILITY:
            return True
        return False
    
    def change_age(self, time_step):
        '''Increase the age of the individual

        As a result it may change its stage (i.e. egg -> larva)

        Args:
            time_step (float)
                Increase in age, in days
        '''
        self.age += time_step
        if self.stage == 'egg' and self.hatching and self.age >= self.time_to_hatch:
            self.stage = 'larva'
        elif self.stage == 'larva' and self.larva and self.age >= self.time_to_pupa:
            self.stage = 'pupa'
        elif self.stage == 'pupa' and self.pupa and self.age >= self.time_to_maturation:
            self.age = 0
            self.stage = 'adult'

    def will_develop(self):
        '''Will this individual develop into an adult?
        
        Returns:
           will_develop (bool)
               Wether the individual will eventually develop into an adult
        '''
        if self.stage == 'egg' and not self.hatching:
            return False
        elif self.stage == 'larva' and not self.larva:
            return False
        else:
            return True

    def is_alive(self):
        '''Is this individual still alive?

        Based on age, stage and development probabilities
        
        Returns:
           is_alive (bool)
               Wether the individual is still alive
        '''
        if self.stage == 'adult' and self.age > self.death:
            return False
        elif self.stage == 'egg' and not self.hatching:
            return False
        elif self.stage == 'larva' and not self.larva:
            return False
        elif self.stage == 'pupa' and not self.pupa:
            return False
        else:
            return True

    def _mendelian(self, genotype):
        if random.random() < 0.5:
            return sorted(genotype)[0]
        else:
            return sorted(genotype)[1]

    def form_gamete1(self, ):
        '''Form a gamete for locus 1 (dsx)
        
        If het. DRIVE homing may happen;
        the genotype at locus 2 may block homing

        Returns:
            genotype (str)
                Genotype at locus 1 for this gamete

        Example:
        >>> i.form_gamete1()
        W
        '''
        if self.hom1:
            return tuple(self.genotype1)[0]
        if DRIVE not in self.genotype1:
            return self._mendelian(self.genotype1)
        else:
            # anti-drive present?
            if ANTI_DRIVE in self.genotype2:
                if random.random() <= ANTI_DRIVE_EFFICIENCY[self.sex]:
                    return self._mendelian(self.genotype1)
            # supermendelian
            if self.nucl_from_father and self.nucl_from_mother:
                if random.random() <= DRIVE_EFFICIENCY_NUCL_FROM_BOTH:
                    return DRIVE
            elif self.nucl_from_father:
                if random.random() <= DRIVE_EFFICIENCY_NUCL_FROM_FATHER.rvs():
                    return DRIVE
            elif self.nucl_from_mother:
                if random.random() <= DRIVE_EFFICIENCY_NUCL_FROM_MOTHER.rvs():
                    return DRIVE
            # resistance
            # TODO: could have functional resistance here too?
            if random.random() < RESISTANCE_EFFICIENCY[self.sex]:
                return RESISTANCE
            # regular mendelian inheritance
            return self._mendelian(self.genotype1)

    def form_gamete2(self):
        '''Form a gamete for locus 2 (antidote)
        
        Returns:
            genotype (str)
                Genotype at locus 2 for this gamete

        Example:
        >>> i.form_gamete2()
        A
        '''
        if self.hom2:
            return tuple(self.genotype2)[0]
        else:
            return self._mendelian(self.genotype2)


def mate_all(population,
             multiple_mating_female=MULTIPLE_MATING_FEMALE,
             multiple_mating_male=MULTIPLE_MATING_MALE):
    '''Randomly mate all adults that can mate

    Args:
        population (iterable)
            All individuals in the population
        multiple_mating_female (bool)
            Wether females can mate multiple times in their lifetime
        multiple_mating_male (bool)
            Wether males can mate multiple times in their lifetime

    Returns:
        eggs (set)
            An iterable of offsprings (Individual objects)
    '''
    eggs = set()
    males = [x for x in population
                if x.sex == 'm'
                and x.stage == 'adult'
                and x.mating
                and ((not multiple_mating_male and not x.mated) or
                     multiple_mating_male)]
    females = [x for x in population
                if x.sex == 'f'
                and x.stage == 'adult'
                and x.mating
                and ((not multiple_mating_female and not x.mated) or
                     multiple_mating_female)]
    random.shuffle(males)
    random.shuffle(females)
    for m, f in zip(males, females):
        eggs = eggs.union(mate(m, f))
    return eggs


def mate(m, f,
         multiple_mating_female=MULTIPLE_MATING_FEMALE,
         multiple_mating_male=MULTIPLE_MATING_MALE):
    '''Mate a female with a male, if conditions are right
    
    Args:
        m (Individual)
            A male
        f (Individual)
            A female
        multiple_mating_female (bool)
            Wether females can mate multiple times in their lifetime
        multiple_mating_male (bool)
            Wether males can mate multiple times in their lifetime

    Returns:
        eggs (set)
            An iterable of offsprings (Individual objects)
    '''
    if f.sex == m.sex:
        raise RuntimeError('Cannot mate')

    eggs = set()

    if not multiple_mating_female:
        f.mated = True
    if not multiple_mating_male:
        m.mated = True

    # can they both mate?
    # also can the female depose eggs?
    if not f.mating or not m.mating or not f.deposing_eggs:
        return eggs

    # we assume all female gametes become eggs
    for egg in range(f.eggs):
        # are we inheriting nuclease from one of the parents?
        nucl_from_father = False
        nucl_from_mother = False
        if DRIVE in f.genotype1:
            nucl_from_mother = True
        if DRIVE in m.genotype1:
            nucl_from_father = True
        # maternal gamete (locus1)
        fg1 = f.form_gamete1()
        # paternal gamete (locus1)
        mg1 = m.form_gamete1()
        # maternal gamete (locus2)
        fg2 = f.form_gamete2()
        # paternal gamete (locus2)
        mg2 = m.form_gamete2()
        # sex
        if random.random() < 0.5:
            sex = 'm'
        else:
            sex = 'f'
        eggs.add(Individual(sex, (fg1, mg1), (fg2, mg2),
                            nucl_from_father, nucl_from_mother))

    # regenerate mating probability for next cycle
    # m.mating = m.get_mating()
    # f.mating = f.get_mating()

    return eggs


def get_all_genotypes():
    '''A generator of all possible genotypes
    
    Useful for output generation'''
    g1 = {DRIVE, RESISTANCE, WILD_TYPE}
    g2 = {ANTI_DRIVE, WILD_TYPE}
    gt1 = set()
    for g11, g12 in itertools.product(g1, g1):
        gt1.add(tuple(sorted((g11, g12))))
    gt2 = set()
    for g21, g22 in itertools.product(g2, g2):
        gt2.add(tuple(sorted((g21, g22))))
    for g1 in sorted(gt1):
        for g2 in sorted(gt2):
            yield ''.join(g1) + ''.join(g2)


def print_header():
    '''Print to stdout the header for the simulation output'''
    print('\t'.join(['round', 'time', 'initial_release',
                     'pop', 'fpop', 'eggs', 'feggs', 'output', 'foutput',
                     'fitness'] +
                    ['WT', 'transgenes', 'drives', 'anti', 'resistance'] +
                    [gt for gt in get_all_genotypes()] +
                    ['%s.eggs' % gt for gt in get_all_genotypes()] +
                    ['WT.output', 'transgenes.output',
                     'drives.output', 'anti.output',
                     'resistance.output'] +
                    ['%s.output' % gt for gt in get_all_genotypes()] 
                    ))


def print_status(time, population, output,
                 initial_population,
                 eggs, repetition):
    '''Print information about the genotype frequencies to stdout
    
    Args:
        time (float)
            Simulation time, in days
        population (iterable)
            All individuals in the population
        output (iterable)
            All larvae + pupae currently available
        initial_population (bool)
            Wether we are introducing the start population
        eggs (iterable)
            All eggs produced at this time point
        repetition (int)
            Round of simulation
    '''
    fpop = [x for x in population if x.sex == 'f']
    pop = len(population)
    results = [repetition, time, initial_population,
               len(population), len(fpop),
               len(eggs), len([x for x in eggs if x.sex == 'f']),
               len(output), len([x for x in output if x.sex == 'f'])]
    if pop != 0:
        if len(fpop) == 0:
            results.append(np.nan)
        else:
            results.append(len([x for x in fpop
                                if x.mating and x.deposing_eggs]) / len(fpop))
        results.append(len([x for x in population
                            if x.genotype1 == {WILD_TYPE, }
                            and x.genotype2 == {WILD_TYPE, }]) / pop)
        results.append(len([x for x in population
                            if x.genotype1.difference({WILD_TYPE, }) != set()
                            or x.genotype2.difference({WILD_TYPE, }) != set()]) / pop)
        results.append(len([x for x in population
                            if DRIVE in x.genotype1]) / pop)
        results.append(len([x for x in population
                            if ANTI_DRIVE in x.genotype2]) / pop)
        results.append(len([x for x in population
                            if RESISTANCE in x.genotype1]) / pop)
        for gt in get_all_genotypes():
            results.append(len([x for x in population
                                if x.get_genotype() == gt]) / pop)
    else:
        for _ in range(6):
            results.append(np.nan)
        for gt in get_all_genotypes():
            results.append(np.nan)
    if len(eggs) != 0:
        for gt in get_all_genotypes():
            results.append(len([x for x in eggs
                                if x.get_genotype() == gt]) / len(eggs))
    else:
        for gt in get_all_genotypes():
            results.append(np.nan)
    if len(output) != 0:
        results.append(len([x for x in output
                            if x.genotype1 == {WILD_TYPE, }
                            and x.genotype2 == {WILD_TYPE, }]) / len(output))
        results.append(len([x for x in output
                            if x.genotype1.difference({WILD_TYPE, }) != set()
                            or x.genotype2.difference({WILD_TYPE, }) != set()]) / len(output))
        results.append(len([x for x in output
                            if DRIVE in x.genotype1]) / len(output))
        results.append(len([x for x in output
                            if ANTI_DRIVE in x.genotype2]) / len(output))
        results.append(len([x for x in output
                            if RESISTANCE in x.genotype1]) / len(output))
        for gt in get_all_genotypes():
            results.append(len([x for x in output
                                if x.get_genotype() == gt]) / len(output))
    else:
        results.append(np.nan)
        results.append(np.nan)
        results.append(np.nan)
        results.append(np.nan)
        results.append(np.nan)
        for gt in get_all_genotypes():
            results.append(np.nan)

    print('\t'.join([str(x) for x in results]))


def run_simulation(start_populations,
                   repetition=0, end_time=365,
                   time_step=TIME_STEP, release=RELEASE,
                   special_releases=None,
                   report_times=None, release_days=RELEASE_DAYS,
                   additional_releases=None,
                   eggs_filter=None):
    '''Run a full large-cage simulation given a series of start populations
    
    Args:
        start_populations (iterable of iterables)
            An iterable of default populations to introduce
            first, alongside the offspring of the whole cage.
            The first element of the iterable is the content of
            the large cage at time point zero.
        repetition (int)
            Round of simulation (useful for reporting)
        end_time (float)
            Maximum length of the simulation (days)
        time_step (float)
            Increase in time each time the simulation moves forward
            NOTE: A time step larger than 0.1 may cause unforeseen bugs
        release (int)
            Maximum number of pupae released on release days
        special_releases (dict)
            key: time, value: release size for that time
        report_times (iterable of int)
            Days for which to report genotype frequencies; by default
            it is done every day
        release_days (iterable of int)
            Days of the week for releases and blood meals. Zero corresponds
            to Monday, six to Sunday
        additional_releases (tuple)
            Additional releases; first element is an iterable of adults,
            the second element is the timepoint at which to start the
            release(s), the third one is the number of releases,
            -1 indicates indefinite releases
        eggs_filter (tuple)
            Trim the eggs output, according to a desired normal distribution
            First element is the loc parameter, second is the scale.
    '''
    if report_times is None:
        report_times = []
    if special_releases is None:
        special_releases = {}
    if eggs_filter is not None:
        eggs_filter_norm = stats.norm(loc=eggs_filter[0],
                                      scale=eggs_filter[1])

    total_time = -time_step

    eggs_nursery = set()
    eggs_hatched = False

    wt_triggered = False

    latest_eggs = set()
    # eggs are harvested in the next feeding cycle
    previous_eggs = set()

    # reverse the order of initial populations
    # so that we can use the "pop" function
    start_populations = start_populations[::-1]
    population = start_populations.pop()

    # counter for additional releases
    additional_releases_counter = 0

    while len(population) > 0 and total_time < end_time:
        initial_population = False
        total_time += time_step
        total_time = round(total_time, 1)
        if total_time == 0:
            initial_population = True

        # age
        dead = set()
        for i in population:
            i.change_age(time_step)
            if not i.is_alive():
                dead.add(i)
        for i in dead:
            population.remove(i)
        # also in egg nursery
        dead = set()
        for e in eggs_nursery:
            e.change_age(time_step)
            if not e.is_alive():
                dead.add(e)
        for e in dead:
            eggs_nursery.remove(e)
        # also in previous egg batch
        dead = set()
        for e in previous_eggs:
            e.change_age(time_step)
            if not e.is_alive():
                dead.add(e)
        for e in dead:
            previous_eggs.remove(e)

        # day of the week
        day = round(total_time, 1) % 7
        
        eggs = set()

        #  Select the larvae from previous harvests
        if len(eggs_nursery) > 0:
            larvae = [x for x in eggs_nursery if x.stage == 'larva']
            pupae = [x for x in eggs_nursery if x.stage == 'pupa']
        else:
            larvae = []
            pupae = []
        
        # feeding/harvesting/release day
        if day % 1 == 0 and int(day) in release_days:
            # collect the previous round of eggs
            if len(previous_eggs) > 0:
                eggs_nursery = eggs_nursery.union(previous_eggs)
            previous_eggs = set()
            latest_eggs = set()

            # add further start populations
            if len(start_populations) > 0 and total_time > 1:
                initial_population = True
                for indv in start_populations.pop():
                    population.add(indv)

            # additional releases (to be done before mating)
            if additional_releases is not None and total_time >= additional_releases[1]:
                additional_releases_counter += 1
                if additional_releases[2] == -1 or additional_releases_counter <= additional_releases[2]:
                    for adult in additional_releases[0]:
                        population.add(deepcopy(adult))

            # mate adults (we are after feeding)
            eggs = mate_all(population)
            # trim eggs if parameter is set
            if eggs_filter is not None:
                eggs = list(eggs)
                random.shuffle(eggs)
                eggs_to_keep = int(eggs_filter_norm.rvs())
                if eggs_to_keep < 0:
                    eggs_to_keep = 0
                elif eggs_to_keep > len(eggs):
                    eggs_to_keep = len(eggs)
                eggs = set(eggs[:eggs_to_keep])

            # save current egg status
            latest_eggs = latest_eggs.union(deepcopy(eggs))
            previous_eggs = previous_eggs.union(deepcopy(eggs))

            if len(pupae) > 0:
                # pick 400 random new pupae to introduce
                random.shuffle(pupae)
                if not round(total_time, 2) % 1 and int(total_time) in special_releases:
                    special = special_releases[int(total_time)]
                    release_pupae = pupae[:special]
                else:
                    release_pupae = pupae[:release]
                population = population.union(release_pupae)
                # remove eggs from nursery
                for p in pupae:
                    eggs_nursery.remove(p)
                eggs_hatched = True

        if ((len(report_times) == 0 and not round(total_time, 2) % 1) or
            round(total_time, 2) in report_times):
            print_status(total_time, population, larvae + pupae,
                         initial_population,
                         latest_eggs, repetition)
