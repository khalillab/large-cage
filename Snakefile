import numpy as np

scenarios = ['base', 'alt', 'lowfitness']
sizes = ['large', 'bugdorm']
cages = ['antidote', 'baseline', 'wt']

rule all:
  input:
    expand('out/{size}/{scenario}/{cage}.tsv',
           scenario=scenarios, cage=cages, size=sizes)

rule simulate:
  input:
      'parameters/{size}/{scenario}/{cage}.yaml'
  output:
      'out/{size}/{scenario}/{cage}.tsv'
  shell:
      '''
      python3 src/simulation.py \
          --parameters {input} \
          > {output}
      '''
