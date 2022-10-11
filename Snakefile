import numpy as np

# specific scenarios
lowfitness = ['lowfitness'] + ['lowfitness_%.2f' % p
                               for p in [0.01, 0.05, 0.2, 0.3, 0.5,
                                         0.7, 0.9, 1]]
releases = [600, 1200, 2400]
times = [[1, 4], [1, 4, 7],
         [0, 1, 2, 3, 4, 5, 6]]
rd = []
for r in releases:
    for t in times:
        if r == 600 and len(t) == 2:
            continue
        rd.append(f'release_{r}_{len(t)}')

# all other scenarios
sizes = ['large', 'bugdorm']
transgenes = ['antidote', 'baseline']
cages = transgenes + ['wt']

rule all:
  input:
    expand('out/{size}/base/{cage}.xlsx',
           cage=cages, size=sizes),
    expand('out/{size}/alt/{cage}.xlsx',
           cage=transgenes, size=sizes),
    expand('out/{size}/{scenario}/antidote.xlsx',
           scenario=lowfitness, size=sizes),
    expand('out/large/{scenario}/{cage}.xlsx',
           scenario=rd, cage=cages)

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

rule convert:
  input:
      'out/{size}/{scenario}/{cage}.tsv'
  output:
      'out/{size}/{scenario}/{cage}.xlsx'
  shell:
      'python3 src/utils/tsv2excel.py {input} {output}'
