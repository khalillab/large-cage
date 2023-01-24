import numpy as np

# specific scenarios
lowfitness = ['lowfitness_%.2f' % p
              for p in [0.01, 0.05, 0.2, 0.3, 0.5,
                        0.7, 0.9, 1]]
releases = range(500, 5100, 500)
drive_thresholds = np.linspace(0.1, 0.8, 15)
drive = []
anti = []
for r in releases:
    drive.append(f'release_{r}')
    for t in drive_thresholds:
        anti.append('release_%d_%.2f' % (r, t))

# all other scenarios
sizes = ['large', 'bugdorm']
antidote = ['antidote']
baseline = ['baseline']
transgenes = antidote + baseline
cages = transgenes + ['wt']

# runs
repeats = 50

rule all:
  input:
    expand('out/{size}/base/{cage}.xlsx',
           cage=cages, size=sizes),
    expand('out/{size}/alt/{cage}.xlsx',
           cage=transgenes, size=sizes),
    expand('out/{size}/{scenario}/antidote.xlsx',
           scenario=lowfitness, size=sizes),
    expand('out/{size}/{scenario}/baseline.xlsx',
           scenario=drive, size=sizes),
    expand('out/{size}/{scenario}/antidote.xlsx',
           scenario=anti, size=sizes)

rule collated:
  input:
    expand('out/{size}/base/{cage}.tsv',
           cage=cages, size=sizes),
    expand('out/{size}/alt/{cage}.tsv',
           cage=transgenes, size=sizes),
    expand('out/{size}/{scenario}/antidote.tsv',
           scenario=lowfitness, size=sizes),
    expand('out/{size}/{scenario}/baseline.tsv',
           scenario=drive, size=sizes),
    expand('out/{size}/{scenario}/antidote.tsv',
           scenario=anti, size=sizes)

rule raw:
  input:
    expand('raw/{size}/base/{cage}_{repeat}.tsv',
           cage=cages, size=sizes, repeat=range(repeats)),
    expand('raw/{size}/alt/{cage}_{repeat}.tsv',
           cage=transgenes, size=sizes, repeat=range(repeats)),
    expand('raw/{size}/{scenario}/antidote_{repeat}.tsv',
           scenario=lowfitness, size=sizes, repeat=range(repeats)),
    expand('raw/{size}/{scenario}/baseline_{repeat}.tsv',
           scenario=drive, size=sizes, repeat=range(repeats)),
    expand('raw/{size}/{scenario}/antidote_{repeat}.xlsx',
           scenario=anti, size=sizes, repeat=range(repeats))

rule simulate:
  input:
      'parameters/{size}/{scenario}/{cage}.yaml'
  output:
      'out/{size}/{scenario}/{cage}_{repeat}.tsv'
  shell:
      '''
      python3 src/simulation.py \
          --parameters {input} \
          > {output}
      '''

rule collate:
  input:
      expand('out/{{size}}/{{scenario}}/{{cage}}_{repeat}.tsv',
             repeat=range(repeats))
  output:
      'out/{size}/{scenario}/{cage}.tsv'
  shell:
      'python3 src/utils/combine_runs.py {input} {output}'

rule convert:
  input:
      'out/{size}/{scenario}/{cage}.tsv'
  output:
      'out/{size}/{scenario}/{cage}.xlsx'
  shell:
      'python3 src/utils/tsv2excel.py {input} {output}'
