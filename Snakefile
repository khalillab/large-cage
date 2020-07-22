POP_SIZE = 400
REPETITIONS = 100
DRIVES = [0.25, 0.5]
ANTIS = [0, 0.25, 0.5, 0.75, 1]
ANTI_EFFECTS = [0.3, 0.6, 0.9]

rule files:
  params:
    time_points = 'data/reference_time_points.txt'

files = rules.files.params

rule baseline:
  message:
    '''
    Run simulations with WT population
    '''
  input:
      time_points = files.time_points
  output:
      'out/wt.tsv'
  params:
      pop_size = POP_SIZE,
      repetitions = REPETITIONS
  shell:
      '''
      python3 src/simulation.py \
          --drive 0 0 \
          --antidote 0 0 \
          --pop-size {params.pop_size} \
          --repetitions {params.repetitions} \
          > {output}
      '''

rule no_antidote:
  message:
    '''
    Run simulations without the antidote

    (custom run based on prelim data)
    '''
  input:
      time_points = files.time_points
  output:
      drive_25 = 'out/no_antidote/drive_25.tsv',
      drive_50 = 'out/no_antidote/drive_50.tsv'
  params:
      pop_size = POP_SIZE,
      repetitions = REPETITIONS,
      drive_25_1 = 71,
      drive_25_2 = 72,
      drive_50_1 = 142,
      drive_50_2 = 143
  shell:
      '''
      python3 src/simulation.py \
          --drive {params.drive_25_1} {params.drive_25_2} \
          --antidote 0 0 \
          --pop-size {params.pop_size} \
          --repetitions {params.repetitions} \
          --time-points {input.time_points} > {output.drive_25}
      python3 src/simulation.py \
          --drive {params.drive_50_1} {params.drive_50_2} \
          --antidote 0 0 \
          --pop-size {params.pop_size} \
          --repetitions {params.repetitions} \
          --time-points {input.time_points} > {output.drive_50}
      '''


rule simulations:
  input:
    expand('out/simulation/{effect}-{drive}-{anti}.tsv',
           effect=ANTI_EFFECTS,
           drive=DRIVES,
           anti=ANTIS)

rule run_simulation:
  message:
    '''
    Run simulations
    '''
  output:
      'out/simulation/{effect}-{drive}-{anti}.tsv'
  params:
      pop_size = POP_SIZE,
      repetitions = REPETITIONS,
      drive = lambda wildcards: int(float(wildcards.drive) * POP_SIZE / 2),
      anti = lambda wildcards: int(float(wildcards.anti) * POP_SIZE / 2)
  shell:
      '''
      python3 src/simulation.py \
          --drive {params.drive} {params.drive} \
          --antidote {params.anti} {params.anti} \
          --pop-size {params.pop_size} \
          --het-antidote-effect {wildcards.effect} \
          --repetitions {params.repetitions} > {output}
      '''

rule simulations_no_antidote_cost:
  input:
    expand('out/simulation_antidote_cost/{hom}-{het}-{drive}-{anti}.tsv',
           het=[1], hom=[1],
           drive=DRIVES,
           anti=ANTIS)

rule run_simulation_no_antidote_cost:
  message:
    '''
    Run simulations (antidote has no fitness cost)
    '''
  output:
      'out/simulation_antidote_cost/{hom}-{het}-{drive}-{anti}.tsv'
  params:
      pop_size = POP_SIZE,
      repetitions = REPETITIONS,
      drive = lambda wildcards: int(float(wildcards.drive) * POP_SIZE / 2),
      anti = lambda wildcards: int(float(wildcards.anti) * POP_SIZE / 2)
  shell:
      '''
      python3 src/simulation.py \
          --drive {params.drive} {params.drive} \
          --antidote {params.anti} {params.anti} \
          --pop-size {params.pop_size} \
          --hom-antidote-effect {wildcards.hom} \
          --het-antidote-effect {wildcards.het} \
          --repetitions {params.repetitions} > {output}
      '''

rule all:
  input:
    rules.baseline.output,
    rules.simulations.input,
    rules.simulations_no_antidote_cost.input,
    rules.no_antidote.output.drive_25,
    rules.no_antidote.output.drive_50
