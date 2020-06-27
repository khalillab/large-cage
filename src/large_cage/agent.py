#!/usr/bin/env python

import sys
import random
import itertools
import numpy as np
from scipy import stats
from copy import deepcopy


# parameters 
POPULATION = 800
SUBSEQUENT_POPULATION = int(POPULATION / 2)
POPULATION_DRIVE = 0.5
POPULATION_ANTI_DRIVE = 0.0
REPETITIONS = 10
HOM_ANTIDRIVE_EFFECT = 0.15
HET_ANTIDRIVE_EFFECT = 0.3
INTRODUCTIONS = 6
GD_INTRODUCTIONS = 6
GD_TRIGGER_DRIVE = 0.9
GD_TRIGGER_DAY = 45
MATING_PROBABILITY = 0.5
MULTIPLE_MATING_FEMALE = False
MULTIPLE_MATING_MALE = True
EGG_DEPOSITION_PROBABILITY = 0.2
ANTI_DRIVE_MATING_HET = MATING_PROBABILITY
ANTI_DRIVE_MATING_HOM = MATING_PROBABILITY
SURVIVAL = stats.weibull_min(c=2.2472084592310644,
                             scale=6.213064151445494,
                             loc=0.7275226197070571)
TIME_STEP = 0.1
NON_FUNCTIONAL = [{'D'}, {'D', 'R'}, {'R'}]
WILD_TYPE = 'W'
DRIVE = 'D'
ANTI_DRIVE = 'A'
RESISTANCE = 'R'
DRIVE_EFFICIENCY_NUCL_FROM_FATHER = stats.norm(loc=0.9867,
                                               scale=0.0163)
DRIVE_EFFICIENCY_NUCL_FROM_MOTHER = stats.norm(loc=0.9667,
                                               scale=0.0408)
DRIVE_EFFICIENCY_NUCL_FROM_BOTH = 1.
ANTI_DRIVE_EFFICIENCY = {'f': 1.,
                         'm': 1.}
RESISTANCE_EFFICIENCY = {'f': 0.4685,
                         'm': 0.4685}
EGGS_WT = stats.norm(loc=7476/150,
                     scale=329/150)
EGGS_NUCL_FROM_MOTHER = stats.norm(loc=2429/150,
                                   scale=472/150)
EGGS_NUCL_FROM_FATHER = stats.norm(loc=5215/150,
                                   scale=447/150)
EGGS_NUCL_FROM_BOTH = stats.norm(loc=3023/150,
                                 scale=522/150)
HATCHING_WT = stats.norm(loc=0.8667,
                         scale=0.0046)
HATCHING_NUCL_FROM_MOTHER = stats.norm(loc=0.5313,
                                       scale=0.0168)
HATCHING_NUCL_FROM_FATHER = stats.norm(loc=0.8725,
                                       scale=0.0159)
HATCHING_NUCL_FROM_BOTH = stats.norm(loc=0.5094,
                                     scale=0.0553)
LARVAL_WT = stats.norm(loc=0.0825,
                       scale=0.0214)
LARVAL_NUCL_FROM_MOTHER = stats.norm(loc=0.1019,
                                     scale=0.0168)
LARVAL_NUCL_FROM_FATHER = stats.norm(loc=0.0671,
                                     scale=0.0146)
LARVAL_NUCL_FROM_BOTH = stats.norm(loc=0.0949,
                                   scale=0.0214)
PUPAL_M_WT = stats.norm(loc=0.1837,
                        scale=0.0627)
PUPAL_M_NUCL_FROM_MOTHER = stats.norm(loc=0.0783,
                                      scale=0.0167)
PUPAL_M_NUCL_FROM_FATHER = stats.norm(loc=0.0692,
                                      scale=0.0233)
PUPAL_M_NUCL_FROM_BOTH = stats.norm(loc=0.0796,
                                    scale=0.0182)
PUPAL_F_WT = stats.norm(loc=0.0918,
                        scale=0.0161)
PUPAL_F_NUCL_FROM_MOTHER = stats.norm(loc=0.0503,
                                      scale=0.0121)
PUPAL_F_NUCL_FROM_FATHER = stats.norm(loc=0.0676,
                                      scale=0.0223)
INTERSEX_NUCL_FROM_FATHER = stats.norm(loc=0.0353,
                                       scale=0.0124)
INTERSEX_NUCL_FROM_MOTHER = stats.norm(loc=0.0096,
                                       scale=0.0066)

class Individual():
    def __init__(self, sex, genotype1, genotype2, x, y,
                 nucl_from_father=False, nucl_from_mother=False):
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
        # around 2 days to hatch
        self.time_to_hatch = random.uniform(1.8, 2.2)
        # around 10 days to pupa
        self.time_to_pupa = random.uniform(9, 10)
        # 10-12 days to mature
        self.time_to_maturation = random.uniform(10, 12)
        # lifespan after full maturation
        self.death = SURVIVAL.rvs()
        # egg -> larva -> pupa -> adult
        self.stage = 'egg'

        # keep track of mating events
        self.mated = False

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

        # position
        self.x = x
        self.y = y

    def get_genotype(self):
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
        if self.sex == 'f' and ANTI_DRIVE in self.genotype2 and self.genotype1 not in NON_FUNCTIONAL:
            if self.hom2:
                return ANTI_DRIVE_MATING_HOM
            else:
                return ANTI_DRIVE_MATING_HET
        else:
            return MATING_PROBABILITY

    def get_mating(self):
        if self._is_intersex():
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
        self.age += time_step
        if self.stage == 'egg' and self.hatching and self.age >= self.time_to_hatch:
            self.stage = 'larva'
        elif self.stage == 'larva' and self.larva and self.age >= self.time_to_pupa:
            self.stage = 'pupa'
        elif self.stage == 'pupa' and self.pupa and self.age >= self.time_to_maturation:
            self.age = 0
            self.stage = 'adult'

    def will_develop(self):
        if self.stage == 'egg' and not self.hatching:
            return False
        elif self.stage == 'larva' and not self.larva:
            return False
        else:
            return True

    def is_alive(self):
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
            return self._mendelian(self.genotype1)

    def form_gamete2(self):
        if self.hom2:
            return tuple(self.genotype2)[0]
        else:
            return self._mendelian(self.genotype2)


def mate_all(population,
             multiple_mating_female=MULTIPLE_MATING_FEMALE,
             multiple_mating_male=MULTIPLE_MATING_MALE):
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
    '''Mate a female with a male, determine the zygotes genotypes'''
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
        eggs.add(Individual(sex, (fg1, mg1), (fg2, mg2), 0, 0,
                            nucl_from_father, nucl_from_mother))

    # regenerate mating probability for next cycle
    # m.mating = m.get_mating()
    # f.mating = f.get_mating()

    return eggs


def get_all_genotypes():
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
    print('\t'.join(['round', 'time', 'initial_release', 'wt_release',
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
                 wt_population, eggs, repetition):
    fpop = [x for x in population if x.sex == 'f']
    pop = len(population)
    results = [repetition, time, initial_population,
               wt_population, len(population), len(fpop),
               len(eggs), len([x for x in eggs if x.sex == 'f']),
               len(output), len([x for x in output if x.sex == 'f'])]
    if pop != 0:
        results.append(len([x for x in fpop
                            if x.mating]) / pop)
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


def run_simulation(start_populations, wt_populations=None,
                   repetition=0, end_time=365, end_drive=1.1,
                   time_step=TIME_STEP, release=SUBSEQUENT_POPULATION,
                   report_times=None):
    if wt_populations is None:
        wt_populations = []
    if report_times is None:
        report_times = []

    total_time = 0

    eggs_nursery = set()
    eggs_hatched = False

    wt_triggered = False

    latest_eggs = set()
    # eggs are harvested in the next feeding cycle
    previous_eggs = set()

    population = start_populations.pop()

    while len(population) > 0 and total_time < end_time:
        initial_population = False
        wt_population = False
        if len(population) > 0:
            drive = len([x for x in population
                         if DRIVE in x.genotype1]) / len(population)
            if drive >= end_drive:
                break
        total_time += time_step
        total_time = round(total_time, 1)

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

        # TODO: move around

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
        if day < .09 or 2.91 < day < 3.09: # or 3.91 < day < 4.09:
            # TODO: spatial effects
            
            # collect the previous round of eggs
            if len(previous_eggs) > 0:
                eggs_nursery = eggs_nursery.union(previous_eggs)
            previous_eggs = set()

            # mate adults (we are after feeding)
            eggs = mate_all(population)

            # save current egg status
            latest_eggs = deepcopy(eggs)
            previous_eggs = deepcopy(eggs)

            if len(start_populations) > 0 and total_time > 1:
                initial_population = True
                for indv in start_populations.pop():
                    population.add(indv)
            if len(pupae) > 0:
                # pick 400 random new pupae to introduce
                random.shuffle(pupae)
                pupae = pupae[:release]
                population = population.union(pupae)
                # remove eggs from nursery
                for p in pupae:
                    eggs_nursery.remove(p)
                eggs_hatched = True
            # subsequent introductions
            # legacy code
            #if (wt_triggered or 
            #    total_time > GD_TRIGGER_DAY or
            #    len([x for x in population
            #         if x.sex == 'f'
            #         and DRIVE in x.genotype1]) / len(population) > GD_TRIGGER_DRIVE):
            #    wt_triggered = True
            #    if len(wt_populations) > 0:
            #        wt_population = True
            #        for indv in wt_populations.pop():
            #            population.add(indv)


        if ((len(report_times) == 0 and not round(total_time, 2) % 1) or
            round(total_time, 2) in report_times):
            #print(total_time, len(population), len(latest_eggs), len(eggs_nursery))
            print_status(total_time, population, larvae + pupae,
                         initial_population,
                         wt_population,
                         latest_eggs, repetition)

