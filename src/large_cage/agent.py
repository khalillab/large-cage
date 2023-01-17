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
params = {
 }

distributions = {}


def get_rvs(dist, params):
    if dist not in distributions:
        distributions[dist] = {}

    if dist == 'norm':
        loc, scale = params
        if params not in distributions[dist]:
            distributions[dist][params] = stats.norm(loc=loc,
                                                     scale=scale)
    elif dist == 'weibull':
        loc, scale, c = params
        if params not in distributions[dist]:
            distributions[dist][params] = stats.weibull_min(c=c,
                                                            loc=loc,
                                                            scale=scale)
    else:
        raise ValueError(f'{dist} not implemented yet')

    return distributions[dist][params].rvs()


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
                 nucl_from_father=False, nucl_from_mother=False,
                 hatching_mod=1, parameters=None):
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
            hatching_mod (float)
                Modifier for the hatching probability
            parameters (dict)
                parameters dictionary
        '''
        if parameters is None:
            self.p = {}
        else:
            self.p = parameters

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
        self.time_to_hatch = random.uniform(self.p['TIME_TO_HATCH'][0], self.p['TIME_TO_HATCH'][1])
        self.time_to_pupa = random.uniform(self.p['TIME_TO_PUPA'][0], self.p['TIME_TO_PUPA'][1])
        self.time_to_maturation = random.uniform(self.p['TIME_TO_MATURATION'][0], self.p['TIME_TO_MATURATION'][1])
        # lifespan after full maturation
        if self.sex == 'm':
            self.death = get_rvs('weibull', (self.p['SURVIVAL_MALE']['loc'],
                                             self.p['SURVIVAL_MALE']['scale'],
                                             self.p['SURVIVAL_MALE']['c']))
        else:
            self.death = get_rvs('weibull', (self.p['SURVIVAL_FEMALE']['loc'],
                                             self.p['SURVIVAL_FEMALE']['scale'],
                                             self.p['SURVIVAL_FEMALE']['c']))

        # egg -> larva -> pupa -> adult
        self.stage = 'egg'

        # keep track of mating events
        self.mated = False

        # is this individual intersex
        self.intersex = self._is_intersex()

        # initial mating probability 
        self.mating = self.get_mating()

        # will deposit eggs?
        # with own genotype modifier
        # when mating will recompute with interacting genotypes
        # we keep this to keep track of population fitness
        self.deposing_eggs = self.deposes_eggs(self.get_deposition_mod())

        # egg production
        if self.sex == 'f' and ''.join(sorted(self.genotype1)) not in [set(x)
                                                                       for x in self.p['NON_FUNCTIONAL']]:
            if self.p['DRIVE'] not in self.genotype1:
                self.eggs = int(get_rvs('norm', (self.p['EGGS_WT']['loc'],
                                                 self.p['EGGS_WT']['scale'])))
            elif self.p['DRIVE'] in self.genotype1 and nucl_from_father and nucl_from_mother:
                self.eggs = int(get_rvs('norm', (self.p['EGGS_NUCL_FROM_BOTH']['loc'],
                                                 self.p['EGGS_NUCL_FROM_BOTH']['scale'])))
            elif self.p['DRIVE'] in self.genotype1 and nucl_from_father:
                self.eggs = int(get_rvs('norm', (self.p['EGGS_NUCL_FROM_FATHER']['loc'],
                                                 self.p['EGGS_NUCL_FROM_FATHER']['scale'])))
            elif self.p['DRIVE'] in self.genotype1 and nucl_from_mother:
                self.eggs = int(get_rvs('norm', (self.p['EGGS_NUCL_FROM_MOTHER']['loc'],
                                                 self.p['EGGS_NUCL_FROM_MOTHER']['scale'])))
            else:
                raise RuntimeError(self.genotype1, nucl_from_father, nucl_from_mother)
        else:
            self.eggs = 0

        if self.p['ANTI_DRIVE'] in self.genotype2:
            if self.hom2:
                self.eggs = int(round(self.eggs * self.p['HOM_ANTIDRIVE_EFFECT']))
            else:
                self.eggs = int(round(self.eggs * self.p['HET_ANTIDRIVE_EFFECT']))

        # hatching probability
        if self.p['DRIVE'] not in self.genotype1:
            hatching = get_rvs('norm', (self.p['HATCHING_WT']['loc'],
                                        self.p['HATCHING_WT']['scale']))
        elif self.p['DRIVE'] in self.genotype1 and nucl_from_father and nucl_from_mother:
            hatching = get_rvs('norm', (self.p['HATCHING_NUCL_FROM_BOTH']['loc'],
                                        self.p['HATCHING_NUCL_FROM_BOTH']['scale']))
        elif self.p['DRIVE'] in self.genotype1 and nucl_from_father:
            hatching = get_rvs('norm', (self.p['HATCHING_NUCL_FROM_FATHER']['loc'],
                                        self.p['HATCHING_NUCL_FROM_FATHER']['scale']))
        elif self.p['DRIVE'] in self.genotype1 and nucl_from_mother:
            hatching = get_rvs('norm', (self.p['HATCHING_NUCL_FROM_MOTHER']['loc'],
                                        self.p['HATCHING_NUCL_FROM_MOTHER']['scale']))
        else:
            raise RuntimeError(self.genotype1, nucl_from_father, nucl_from_mother)
        # apply modifier for hatching probability
        hatching = hatching * hatching_mod
        if random.random() < hatching:
            self.hatching = True
        else:
            self.hatching = False

        # larval mortality
        if self.p['DRIVE'] not in self.genotype1:
            larval = get_rvs('norm', (self.p['LARVAL_WT']['loc'],
                                      self.p['LARVAL_WT']['scale']))
        elif self.p['DRIVE'] in self.genotype1 and nucl_from_father and nucl_from_mother:
            larval = get_rvs('norm', (self.p['LARVAL_NUCL_FROM_BOTH']['loc'],
                                      self.p['LARVAL_NUCL_FROM_BOTH']['scale']))
        elif self.p['DRIVE'] in self.genotype1 and nucl_from_father:
            larval = get_rvs('norm', (self.p['LARVAL_NUCL_FROM_FATHER']['loc'],
                                      self.p['LARVAL_NUCL_FROM_FATHER']['scale']))
        elif self.p['DRIVE'] in self.genotype1 and nucl_from_mother:
            larval = get_rvs('norm', (self.p['LARVAL_NUCL_FROM_MOTHER']['loc'],
                                      self.p['LARVAL_NUCL_FROM_MOTHER']['scale']))
        else:
            raise RuntimeError(self.genotype1, nucl_from_father, nucl_from_mother)
        if random.random() < larval:
            self.larva = False
        else:
            self.larva = True

        # pupal mortality
        if sex == 'm':
            if self.p['DRIVE'] not in self.genotype1:
                pupal = get_rvs('norm', (self.p['PUPAL_M_WT']['loc'],
                                         self.p['PUPAL_M_WT']['scale']))
            elif self.p['DRIVE'] in self.genotype1 and nucl_from_father and nucl_from_mother:
                pupal = get_rvs('norm', (self.p['PUPAL_M_NUCL_FROM_BOTH']['loc'],
                                         self.p['PUPAL_M_NUCL_FROM_BOTH']['scale']))
            elif self.p['DRIVE'] in self.genotype1 and nucl_from_father:
                pupal = get_rvs('norm', (self.p['PUPAL_M_NUCL_FROM_FATHER']['loc'],
                                         self.p['PUPAL_M_NUCL_FROM_FATHER']['scale']))
            elif self.p['DRIVE'] in self.genotype1 and nucl_from_mother:
                pupal = get_rvs('norm', (self.p['PUPAL_M_NUCL_FROM_MOTHER']['loc'],
                                         self.p['PUPAL_M_NUCL_FROM_MOTHER']['scale']))
            else:
                raise RuntimeError(self.genotype1, nucl_from_father, nucl_from_mother)
        else:
            if self.p['DRIVE'] not in self.genotype1:
                pupal = get_rvs('norm', (self.p['PUPAL_F_WT']['loc'],
                                         self.p['PUPAL_F_WT']['scale']))
            elif self.p['DRIVE'] in self.genotype1 and nucl_from_father and nucl_from_mother:
                # assumed
                pupal = get_rvs('norm', (self.p['PUPAL_M_NUCL_FROM_BOTH']['loc'],
                                         self.p['PUPAL_M_NUCL_FROM_BOTH']['scale']))
            elif self.p['DRIVE'] in self.genotype1 and nucl_from_father:
                pupal = get_rvs('norm', (self.p['PUPAL_F_NUCL_FROM_FATHER']['loc'],
                                         self.p['PUPAL_F_NUCL_FROM_FATHER']['scale']))
            elif self.p['DRIVE'] in self.genotype1 and nucl_from_mother:
                pupal = get_rvs('norm', (self.p['PUPAL_F_NUCL_FROM_MOTHER']['loc'],
                                         self.p['PUPAL_F_NUCL_FROM_MOTHER']['scale']))
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
        if self.p['DRIVE'] in self.genotype1 and (self.nucl_from_father and self.nucl_from_mother):
            intersex = 0
        elif self.p['DRIVE'] in self.genotype1 and self.nucl_from_father and self.sex == 'f':
            intersex = get_rvs('norm', (self.p['INTERSEX_NUCL_FROM_FATHER']['loc'],
                                        self.p['INTERSEX_NUCL_FROM_FATHER']['scale']))
        elif self.p['DRIVE'] in self.genotype1 and self.nucl_from_mother and self.sex == 'f':
            intersex = get_rvs('norm', (self.p['INTERSEX_NUCL_FROM_MOTHER']['loc'],
                                        self.p['INTERSEX_NUCL_FROM_MOTHER']['scale']))
        else:
            intersex = 0
        if random.random() <= intersex:
            return True
        return False

    def _mating_probability(self):
        if self.sex == 'm':
            mating_probability = self.p['MATING_PROBABILITY_MALE']
        else:
            mating_probability = self.p['MATING_PROBABILITY']
        return mating_probability * self.p['MATING_MOD'][self.sex
                ].get(self.get_genotype(), 1)

    def get_mating(self):
        '''Generate the probability that this individual will mate

        Returns
            mating (bool)
               Wether the individual will mate

        The decision is based on sex, genotype and intersex phenotype
        '''
        if self.intersex:
            return False
        if self.sex == 'f' and self.genotype1 in [set(x) for x in
                                                  self.p['NON_FUNCTIONAL']]:
            return False
        else:
            prob = self._mating_probability()
            if random.random() < prob:
                return True
            else:
                return False

    def deposes_eggs(self, modifier):
        '''Generate the probability that this individual will depose eggs

        Args:
            modifier (float)
               Modifier for the eggs deposition modifier

        Returns
            deposes (bool)
               Wether the individual will depose eggs

        The decision is based on sex and modifier
        '''
        if self.sex == 'm':
            return None
        elif random.random() <= self.p['EGG_DEPOSITION_PROBABILITY'] * modifier:
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
        if self.p['DRIVE'] not in self.genotype1:
            return self._mendelian(self.genotype1)
        else:
            # anti-drive present?
            if self.p['ANTI_DRIVE'] in self.genotype2:
                if random.random() <= 1 - self.p['DRIVE_EFFICIENCY_MOD'][self.sex
                        ].get(self.get_genotype(), 1):
                    return self.p['DRIVE']
                else:
                    # return the other allele
                    try:
                        return sorted(set(self.genotype1).difference((self.p['DRIVE'],)))[0]
                    except Exception as e:
                        raise RuntimeError(self.genotype1, self.genotype2, str(e))
            # supermendelian
            if self.sex == 'm':
                if random.random() <= get_rvs('norm', (self.p['DRIVE_EFFICIENCY_MALE']['loc'],
                                                       self.p['DRIVE_EFFICIENCY_MALE']['scale'])):
                    return self.p['DRIVE']
            else:
                if random.random() <= get_rvs('norm', (self.p['DRIVE_EFFICIENCY_FEMALE']['loc'],
                                                       self.p['DRIVE_EFFICIENCY_FEMALE']['scale'])):
                    return self.p['DRIVE']
            # resistance
            # TODO: could have functional resistance here too?
            if random.random() < self.p['RESISTANCE_EFFICIENCY'][self.sex]:
                return self.p['RESISTANCE']
            # return the other allele
            try:
                return sorted(set(self.genotype1).difference((self.p['DRIVE'])))[0]
            except Exception as e:
                raise RuntimeError(self.genotype1, self.genotype2, str(e))

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
            # anti-drive present?
            if self.p['ANTI_DRIVE'] in self.genotype2:
                if random.random() <= self.p['ANTIDOTE_INHERITANCE'][self.sex
                        ].get(self.get_genotype(), 1):
                    return self.p['ANTI_DRIVE']
                else:
                    # return the other allele
                    try:
                        return sorted(set(self.genotype1).difference((self.p['ANTI_DRIVE'],)))[0]
                    except Exception as e:
                        raise RuntimeError(self.genotype1, self.genotype2, str(e))
            # simple mendelian
            return self._mendelian(self.genotype2)

    def get_egg_mod(self):
        '''Get the genotype-specific modifier for the number of eggs

        Returns:
            modifier (float)
                Modifier for this genotype

        Example:
        >>> i.get_egg_mod()
        0.95
        '''
        return self.p['EGGS_MOD'][self.sex].get(self.get_genotype(), 1.)

    def get_hatching_mod(self):
        '''Get the genotype-specific modifier for the hatching rate

        Returns:
            modifier (float)
                Modifier for this genotype

        Example:
        >>> i.get_hatching_mod()
        0.95
        '''
        return self.p['HATCHING_MOD'][self.sex].get(self.get_genotype(), 1.)

    def get_deposition_mod(self):
        '''Get the genotype-specific modifier for the deposition prob

        Returns:
            modifier (float)
                Modifier for this genotype

        Example:
        >>> i.get_deposition_mod()
        0.95
        '''
        return self.p['DEPOSITION_MOD'][self.sex].get(self.get_genotype(), 1.)


def mate_all(population, p=None,
             multiple_mating_female=None,
             multiple_mating_male=None):
    '''Randomly mate all adults that can mate

    Args:
        population (iterable)
            All individuals in the population
        p (dict)
            Parameters
        multiple_mating_female (bool)
            Wether females can mate multiple times in their lifetime
        multiple_mating_male (bool)
            Wether males can mate multiple times in their lifetime

    Returns:
        eggs (set)
            An iterable of offsprings (Individual objects)
    '''
    if p is None:
        p = params
    if multiple_mating_female is None:
        multiple_mating_female=p['MULTIPLE_MATING_FEMALE']
    if multiple_mating_male is None:
        multiple_mating_male=p['MULTIPLE_MATING_MALE']
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
        eggs = eggs.union(mate(m, f, p=p))
    return eggs


def mate(m, f, p=None,
         multiple_mating_female=None,
         multiple_mating_male=None):
    '''Mate a female with a male, if conditions are right

    Args:
        m (Individual)
            A male
        f (Individual)
            A female
        p (dict)
            Parameters
        multiple_mating_female (bool)
            Wether females can mate multiple times in their lifetime
        multiple_mating_male (bool)
            Wether males can mate multiple times in their lifetime

    Returns:
        eggs (set)
            An iterable of offsprings (Individual objects)
    '''
    if p is None:
        p = params
    if multiple_mating_female is None:
        multiple_mating_female=p['MULTIPLE_MATING_FEMALE']
    if multiple_mating_male is None:
        multiple_mating_male=p['MULTIPLE_MATING_MALE']
    if f.sex == m.sex:
        raise RuntimeError('Cannot mate')

    eggs = set()

    if not multiple_mating_female:
        f.mated = True
    if not multiple_mating_male:
        m.mated = True

    # introduce modifiers
    # needs to be done here because male and female genotypes interact
    deposition_mod = m.get_deposition_mod() * f.get_deposition_mod()
    egg_mod = m.get_egg_mod() * f.get_egg_mod()
    hatching_mod = m.get_hatching_mod() * f.get_hatching_mod()

    # can they both mate?
    # also can the female depose eggs?
    if not f.mating or not m.mating:
        return eggs
    if not f.deposes_eggs(deposition_mod):
        return eggs

    actual_eggs = round(f.eggs * egg_mod)

    # we assume all female gametes become eggs
    for egg in range(actual_eggs):
        # are we inheriting nuclease from one of the parents?
        nucl_from_father = False
        nucl_from_mother = False
        if p['DRIVE'] in f.genotype1:
            nucl_from_mother = True
        if p['DRIVE'] in m.genotype1:
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
                            nucl_from_father, nucl_from_mother,
                            hatching_mod, parameters=p))

    # regenerate mating probability for next cycle
    # m.mating = m.get_mating()
    # f.mating = f.get_mating()

    return eggs


def get_all_genotypes(p=None):
    '''A generator of all possible genotypes

    Useful for output generation

    Args:
        p (dict)
            Parameters
    '''
    if p is None:
        p = params
    g1 = {p['DRIVE'], p['RESISTANCE'], p['WILD_TYPE']}
    g2 = {p['ANTI_DRIVE'], p['WILD_TYPE']}
    gt1 = set()
    for g11, g12 in itertools.product(g1, g1):
        gt1.add(tuple(sorted((g11, g12))))
    gt2 = set()
    for g21, g22 in itertools.product(g2, g2):
        gt2.add(tuple(sorted((g21, g22))))
    for g1 in sorted(gt1):
        for g2 in sorted(gt2):
            yield ''.join(g1) + ''.join(g2)


def print_header(p=None):
    '''Print to stdout the header for the simulation output

    Args:
        p (dict)
            Parameters
    '''
    if p is None:
        p = params
    print('\t'.join(['round', 'time', 'initial_release',
                     'pop', 'fpop', 'eggs', 'feggs', 'output', 'foutput',
                     'fitness'] +
                    ['WT', 'transgenes', 'drives', 'anti', 'resistance'] +
                    [gt for gt in get_all_genotypes(p)] +
                    ['%s.eggs' % gt for gt in get_all_genotypes(p)] +
                    ['WT.output', 'transgenes.output',
                     'drives.output', 'anti.output',
                     'resistance.output'] +
                    ['%s.output' % gt for gt in get_all_genotypes(p)] 
                    ))


def get_drive_frequency(population, p=None):
    '''Get frequency of drive individuals in the population

    Args:
        population (iterable)
            All individuals in the population
        p (dict)
            Parameters

    Returns:
        proportion (float)
            The proportion of individuals that carry a drive allele
    '''
    if p is None:
        p = params

    pop = len(population)
    if pop != 0:
        prop = len([x for x in population
                    if p['DRIVE'] in x.genotype1]) / pop
    else:
        prop = np.nan

    return prop


def print_status(time, population, output,
                 initial_population,
                 eggs, repetition, p=None):
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
        p (dict)
            Parameters
    '''
    if p is None:
        p = params
    fpop = [x for x in population if x.sex == 'f']
    pop = len(population)
    results = [repetition, time, initial_population,
               len(population), len(fpop),
               len(eggs), len([x for x in eggs if x.sex == 'f']),
               len(output), len([x for x in output if x.sex == 'f'])]
    if pop != 0:
        if len(fpop) == 0:
            results.append(0.)
        else:
            results.append(len([x for x in fpop
                                if x.mating and x.deposing_eggs]) / len(fpop))
        results.append(len([x for x in population
                            if x.genotype1 == {p['WILD_TYPE'], }
                            and x.genotype2 == {p['WILD_TYPE'], }]) / pop)
        results.append(len([x for x in population
                            if x.genotype1.difference({p['WILD_TYPE'], }) != set()
                            or x.genotype2.difference({p['WILD_TYPE'], }) != set()]) / pop)
        results.append(len([x for x in population
                            if p['DRIVE'] in x.genotype1]) / pop)
        results.append(len([x for x in population
                            if p['ANTI_DRIVE'] in x.genotype2]) / pop)
        results.append(len([x for x in population
                            if p['RESISTANCE'] in x.genotype1]) / pop)
        for gt in get_all_genotypes(p):
            results.append(len([x for x in population
                                if x.get_genotype() == gt]) / pop)
    else:
        for _ in range(6):
            results.append(np.nan)
        for gt in get_all_genotypes(p):
            results.append(np.nan)
    if len(eggs) != 0:
        for gt in get_all_genotypes(p):
            results.append(len([x for x in eggs
                                if x.get_genotype() == gt]) / len(eggs))
    else:
        for gt in get_all_genotypes(p):
            results.append(np.nan)
    if len(output) != 0:
        results.append(len([x for x in output
                            if x.genotype1 == {p['WILD_TYPE'], }
                            and x.genotype2 == {p['WILD_TYPE'], }]) / len(output))
        results.append(len([x for x in output
                            if x.genotype1.difference({p['WILD_TYPE'], }) != set()
                            or x.genotype2.difference({p['WILD_TYPE'], }) != set()]) / len(output))
        results.append(len([x for x in output
                            if p['DRIVE'] in x.genotype1]) / len(output))
        results.append(len([x for x in output
                            if p['ANTI_DRIVE'] in x.genotype2]) / len(output))
        results.append(len([x for x in output
                            if p['RESISTANCE'] in x.genotype1]) / len(output))
        for gt in get_all_genotypes(p):
            results.append(len([x for x in output
                                if x.get_genotype() == gt]) / len(output))
    else:
        results.append(np.nan)
        results.append(np.nan)
        results.append(np.nan)
        results.append(np.nan)
        results.append(np.nan)
        for gt in get_all_genotypes(p):
            results.append(np.nan)

    print('\t'.join([str(x) for x in results]))


def run_simulation(start_populations,
                   repetition=0, end_time=365,
                   time_step=None, release=None,
                   special_releases=None,
                   report_times=None, release_days=None,
                   additional_releases=None,
                   eggs_filter=None,
                   use_adults_if_needed=False,
                   p=None):
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
            -1 indicates indefinite releases. The fourth is the frequency
            of drive individuals at which to start the releases, which
            overrides the second argument if not None
        eggs_filter (tuple)
            Trim the eggs output, according to a desired normal distribution
            First element is the loc parameter, second is the scale.
        use_adults_if_needed (bool)
            If there are no pupae in the egg nursery, use adults
            (can be necessary if there is a single release)
        p (dict)
            Parameters
    '''
    if p is None:
        p = params
    if time_step is None:
        time_step=p['TIME_STEP']
    if release is None:
        release=p['RELEASE']
    if release_days is None:
        release_days=p['RELEASE_DAYS']
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

    # keep track of drive frequencies
    drive_frequencies = []
    drive_ever_released = False
    drive_threshold_passed = False

    while len(population) > 0 and total_time < end_time:
        restocking = False
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
            if len(pupae) == 0 and use_adults_if_needed:
                # could happen if there is a single release day
                pupae = [x for x in eggs_nursery if x.stage == 'adult']
        else:
            larvae = []
            pupae = []

        # feeding/harvesting/release day
        if day % 1 == 0 and int(day) in release_days:
            restocking = True
            # collect the previous round of eggs
            if len(previous_eggs) > 0:
                eggs_nursery = eggs_nursery.union(previous_eggs)
            previous_eggs = set()
            latest_eggs = set()

            # add further start populations
            if len(start_populations) > 0 and total_time > 1:
                for indv in start_populations.pop():
                    population.add(indv)

            # additional releases (to be done before mating)
            if additional_releases is not None:
                if (additional_releases[3] is None and total_time >= additional_releases[1]) or (additional_releases[3] is not None and drive_threshold_passed):
                    additional_releases_counter += 1
                    if additional_releases[2] == -1 or additional_releases_counter <= additional_releases[2]:
                        for adult in additional_releases[0]:
                            population.add(deepcopy(adult))

            # mate adults (we are after feeding)
            eggs = mate_all(population, p=p)
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
                for ep in release_pupae:
                    eggs_nursery.remove(ep)
                eggs_hatched = True

        if ((len(report_times) == 0 and not round(total_time, 2) % 1) or
            round(total_time, 2) in report_times):
            if not drive_threshold_passed:
                # keep track of drive frequencies
                drive_freq = get_drive_frequency(population, p)
                if not drive_ever_released and drive_freq > 0:
                    drive_ever_released = True
                    sys.stderr.write(f'{total_time} drive observed\n')
                drive_frequencies.append(drive_freq)
                drive_frequencies = drive_frequencies[-7:]
                # if set, check if drive frequency threshold has been passed
                if additional_releases[3] is not None and drive_ever_released:
                    if len([x for x in drive_frequencies
                            if x > additional_releases[3]]) == len(drive_frequencies):
                                drive_threshold_passed = True
                                sys.stderr.write(f'{total_time} will start antidote releases\n')

            print_status(total_time, population, larvae + pupae,
                         restocking,
                         latest_eggs, repetition, p)
