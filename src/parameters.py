#!/usr/bin/env python

from scipy import stats

# parameters 
#
# NOTE: change values here and call the simulation.py script
#       with the --override-parameters src/parameters.py argument
#
# NOTE: if a parameter is removed from this file, the default in
#       src/large_cage/agent.py will be used
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
# mating probability antidrive
ANTI_DRIVE_MATING_HET = MATING_PROBABILITY
ANTI_DRIVE_MATING_HOM = MATING_PROBABILITY
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
EGGS_WT = stats.norm(loc=7476/150,
                     scale=329/150)
EGGS_NUCL_FROM_MOTHER = stats.norm(loc=2429/150,
                                   scale=472/150)
EGGS_NUCL_FROM_FATHER = stats.norm(loc=5215/150,
                                   scale=447/150)
EGGS_NUCL_FROM_BOTH = stats.norm(loc=3023/150,
                                 scale=522/150)
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
