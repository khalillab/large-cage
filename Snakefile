POP_SIZE = 400
REPETITIONS = 50
DRIVES = [0.25, 0.5]
ANTIS = [0, 0.5, 1]
ANTI_EFFECTS = [0.3, 0.9]
# late releases
LATE_WT = 100
LATE_ANTI = 100
LATE_START = 45
LATE_COUNT = 10

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
      repetitions = int(REPETITIONS / 5)
  shell:
      '''
      python3 src/simulation.py \
          --drive 0 0 0 0 \
          --antidote 0 0 0 0 \
          --wild-type {params.pop_size} {params.pop_size} {params.pop_size} {params.pop_size} \
          --release {params.pop_size} \
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
          --drive 0 0 0 0 {params.drive_25_1} {params.drive_25_2} \
          --antidote 0 0 0 0 0 0 \
          --wild-type {params.pop_size} {params.pop_size} {params.pop_size} {params.pop_size} 0 0 \
          --release {params.pop_size} \
          --repetitions {params.repetitions} \
          --time-points {input.time_points} > {output.drive_25}
      python3 src/simulation.py \
          --drive 0 0 0 0 {params.drive_50_1} {params.drive_50_2} \
          --antidote 0 0 0 0 0 0 \
          --wild-type {params.pop_size} {params.pop_size} {params.pop_size} {params.pop_size} 0 0 \
          --release {params.pop_size} \
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
          --drive 0 0 0 0 {params.drive} {params.drive} \
          --antidote 0 0 0 0 {params.anti} {params.anti} \
          --wild-type {params.pop_size} {params.pop_size} {params.pop_size} {params.pop_size} 0 0 \
          --release {params.pop_size} \
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
          --drive 0 0 0 0 {params.drive} {params.drive} \
          --antidote 0 0 0 0 {params.anti} {params.anti} \
          --wild-type {params.pop_size} {params.pop_size} {params.pop_size} {params.pop_size} 0 0 \
          --release {params.pop_size} \
          --hom-antidote-effect {wildcards.hom} \
          --het-antidote-effect {wildcards.het} \
          --repetitions {params.repetitions} > {output}
      '''

rule simulations_10:
  input:
    expand('out/simulation_10/{effect}-{drive}-{anti}.tsv',
           effect=ANTI_EFFECTS,
           drive=DRIVES,
           anti=ANTIS)

rule run_simulation_10:
  message:
    '''
    Run simulations (10 releases)
    '''
  output:
      'out/simulation_10/{effect}-{drive}-{anti}.tsv'
  params:
      pop_size = POP_SIZE,
      repetitions = REPETITIONS,
      drive = lambda wildcards: int(float(wildcards.drive) * POP_SIZE / 2),
      anti = lambda wildcards: int(float(wildcards.anti) * POP_SIZE / 2)
  shell:
      '''
      python3 src/simulation.py \
          --drive 0 0 0 0 {params.drive} {params.drive} {params.drive} {params.drive} {params.drive} {params.drive} {params.drive} {params.drive} {params.drive} {params.drive}\
          --antidote 0 0 0 0 {params.anti} {params.anti} {params.anti} {params.anti} {params.anti} {params.anti} {params.anti} {params.anti} {params.anti} {params.anti}\
          --wild-type {params.pop_size} {params.pop_size} {params.pop_size} {params.pop_size} 0 0 0 0 0 0 0 0 0 0 \
          --release {params.pop_size} \
          --het-antidote-effect {wildcards.effect} \
          --repetitions {params.repetitions} > {output}
      '''

rule simulations_10_no_antidote_cost:
  input:
    expand('out/simulation_10_antidote_cost/{hom}-{het}-{drive}-{anti}.tsv',
           het=[1], hom=[1],
           drive=DRIVES,
           anti=ANTIS)

rule run_simulation_10_no_antidote_cost:
  message:
    '''
    Run simulations (10 releases, antidote has no fitness cost)
    '''
  output:
      'out/simulation_10_antidote_cost/{hom}-{het}-{drive}-{anti}.tsv'
  params:
      pop_size = POP_SIZE,
      repetitions = REPETITIONS,
      drive = lambda wildcards: int(float(wildcards.drive) * POP_SIZE / 2),
      anti = lambda wildcards: int(float(wildcards.anti) * POP_SIZE / 2)
  shell:
      '''
      python3 src/simulation.py \
          --drive 0 0 0 0 {params.drive} {params.drive} {params.drive} {params.drive} {params.drive} {params.drive} {params.drive} {params.drive} {params.drive} {params.drive} \
          --antidote 0 0 0 0 {params.anti} {params.anti} {params.anti} {params.anti} {params.anti} {params.anti} {params.anti} {params.anti} {params.anti} {params.anti} \
          --wild-type {params.pop_size} {params.pop_size} {params.pop_size} {params.pop_size} 0 0 0 0 0 0 0 0 0 0 \
          --release {params.pop_size} \
          --hom-antidote-effect {wildcards.hom} \
          --het-antidote-effect {wildcards.het} \
          --repetitions {params.repetitions} > {output}
      '''

rule simulations_late_wt:
  input:
    expand('out/simulation_late_wt/{effect}-{drive}-{anti}.tsv',
           effect=ANTI_EFFECTS,
           drive=DRIVES,
           anti=ANTIS)

rule run_simulation_late_wt:
  message:
    '''
    Run simulations (late WT release)
    '''
  output:
      'out/simulation_late_wt/{effect}-{drive}-{anti}.tsv'
  params:
      pop_size = POP_SIZE,
      repetitions = REPETITIONS,
      drive = lambda wildcards: int(float(wildcards.drive) * POP_SIZE / 2),
      anti = lambda wildcards: int(float(wildcards.anti) * POP_SIZE / 2),
      late_start = LATE_START,
      late_wt = LATE_WT
  shell:
      '''
      python3 src/simulation.py \
          --drive 0 0 0 0 {params.drive} {params.drive} \
          --antidote 0 0 0 0 {params.anti} {params.anti} \
          --wild-type {params.pop_size} {params.pop_size} {params.pop_size} {params.pop_size} 0 0 \
          --release {params.pop_size} \
          --het-antidote-effect {wildcards.effect} \
          --late-releases -1 \
          --late-releases-start {params.late_start} \
          --late-wt {params.late_wt} \
          --repetitions {params.repetitions} > {output}
      '''

rule simulations_late_wt_no_antidote_cost:
  input:
    expand('out/simulation_late_wt_antidote_cost/{hom}-{het}-{drive}-{anti}.tsv',
           het=[1], hom=[1],
           drive=DRIVES,
           anti=ANTIS)

rule run_simulation_late_wt_no_antidote_cost:
  message:
    '''
    Run simulations (late WT release, antidote has no fitness cost)
    '''
  output:
      'out/simulation_late_wt_antidote_cost/{hom}-{het}-{drive}-{anti}.tsv'
  params:
      pop_size = POP_SIZE,
      repetitions = REPETITIONS,
      drive = lambda wildcards: int(float(wildcards.drive) * POP_SIZE / 2),
      anti = lambda wildcards: int(float(wildcards.anti) * POP_SIZE / 2),
      late_start = LATE_START,
      late_wt = LATE_WT
  shell:
      '''
      python3 src/simulation.py \
          --drive 0 0 0 0 {params.drive} {params.drive} \
          --antidote 0 0 0 0 {params.anti} {params.anti} \
          --wild-type {params.pop_size} {params.pop_size} {params.pop_size} {params.pop_size} 0 0 \
          --release {params.pop_size} \
          --hom-antidote-effect {wildcards.hom} \
          --het-antidote-effect {wildcards.het} \
          --late-releases -1 \
          --late-releases-start {params.late_start} \
          --late-wt {params.late_wt} \
          --repetitions {params.repetitions} > {output}
      '''

rule simulations_late_anti:
  input:
    expand('out/simulation_late_anti/{effect}-{drive}-{anti}.tsv',
           effect=ANTI_EFFECTS,
           drive=DRIVES,
           anti=ANTIS)

rule run_simulation_late_anti:
  message:
    '''
    Run simulations (late antidote release)
    '''
  output:
      'out/simulation_late_anti/{effect}-{drive}-{anti}.tsv'
  params:
      pop_size = POP_SIZE,
      repetitions = REPETITIONS,
      drive = lambda wildcards: int(float(wildcards.drive) * POP_SIZE / 2),
      anti = lambda wildcards: int(float(wildcards.anti) * POP_SIZE / 2),
      late_start = LATE_START,
      late_anti = LATE_ANTI
  shell:
      '''
      python3 src/simulation.py \
          --drive 0 0 0 0 {params.drive} {params.drive} \
          --antidote 0 0 0 0 {params.anti} {params.anti} \
          --wild-type {params.pop_size} {params.pop_size} {params.pop_size} {params.pop_size} 0 0 \
          --release {params.pop_size} \
          --het-antidote-effect {wildcards.effect} \
          --late-releases -1 \
          --late-releases-start {params.late_start} \
          --late-anti {params.late_anti} \
          --repetitions {params.repetitions} > {output}
      '''

rule simulations_late_anti_no_antidote_cost:
  input:
    expand('out/simulation_late_anti_antidote_cost/{hom}-{het}-{drive}-{anti}.tsv',
           het=[1], hom=[1],
           drive=DRIVES,
           anti=ANTIS)

rule run_simulation_late_anti_no_antidote_cost:
  message:
    '''
    Run simulations (late antidote release, antidote has no fitness cost)
    '''
  output:
      'out/simulation_late_anti_antidote_cost/{hom}-{het}-{drive}-{anti}.tsv'
  params:
      pop_size = POP_SIZE,
      repetitions = REPETITIONS,
      drive = lambda wildcards: int(float(wildcards.drive) * POP_SIZE / 2),
      anti = lambda wildcards: int(float(wildcards.anti) * POP_SIZE / 2),
      late_start = LATE_START,
      late_anti = LATE_ANTI
  shell:
      '''
      python3 src/simulation.py \
          --drive 0 0 0 0 {params.drive} {params.drive} \
          --antidote 0 0 0 0 {params.anti} {params.anti} \
          --wild-type {params.pop_size} {params.pop_size} {params.pop_size} {params.pop_size} 0 0 \
          --release {params.pop_size} \
          --hom-antidote-effect {wildcards.hom} \
          --het-antidote-effect {wildcards.het} \
          --late-releases -1 \
          --late-releases-start {params.late_start} \
          --late-anti {params.late_anti} \
          --repetitions {params.repetitions} > {output}
      '''

rule simulations_late_anti_once:
  input:
    expand('out/simulation_late_anti_once/{effect}-{drive}-{anti}.tsv',
           effect=ANTI_EFFECTS,
           drive=DRIVES,
           anti=ANTIS)

rule run_simulation_late_anti_once:
  message:
    '''
    Run simulations (late antidote release, not continuously)
    '''
  output:
      'out/simulation_late_anti_once/{effect}-{drive}-{anti}.tsv'
  params:
      pop_size = POP_SIZE,
      repetitions = REPETITIONS,
      drive = lambda wildcards: int(float(wildcards.drive) * POP_SIZE / 2),
      anti = lambda wildcards: int(float(wildcards.anti) * POP_SIZE / 2),
      late_start = LATE_START,
      late_anti = LATE_ANTI,
      late_count = LATE_COUNT
  shell:
      '''
      python3 src/simulation.py \
          --drive 0 0 0 0 {params.drive} {params.drive} \
          --antidote 0 0 0 0 {params.anti} {params.anti} \
          --wild-type {params.pop_size} {params.pop_size} {params.pop_size} {params.pop_size} 0 0 \
          --release {params.pop_size} \
          --het-antidote-effect {wildcards.effect} \
          --late-releases {params.late_count} \
          --late-releases-start {params.late_start} \
          --late-anti {params.late_anti} \
          --repetitions {params.repetitions} > {output}
      '''

rule simulations_late_anti_once_no_antidote_cost:
  input:
    expand('out/simulation_late_anti_once_antidote_cost/{hom}-{het}-{drive}-{anti}.tsv',
           het=[1], hom=[1],
           drive=DRIVES,
           anti=ANTIS)

rule run_simulation_late_anti_once_no_antidote_cost:
  message:
    '''
    Run simulations (late antidote release, not continuously, antidote has no fitness cost)
    '''
  output:
      'out/simulation_late_anti_once_antidote_cost/{hom}-{het}-{drive}-{anti}.tsv'
  params:
      pop_size = POP_SIZE,
      repetitions = REPETITIONS,
      drive = lambda wildcards: int(float(wildcards.drive) * POP_SIZE / 2),
      anti = lambda wildcards: int(float(wildcards.anti) * POP_SIZE / 2),
      late_start = LATE_START,
      late_anti = LATE_ANTI,
      late_count = LATE_COUNT
  shell:
      '''
      python3 src/simulation.py \
          --drive 0 0 0 0 {params.drive} {params.drive} \
          --antidote 0 0 0 0 {params.anti} {params.anti} \
          --wild-type {params.pop_size} {params.pop_size} {params.pop_size} {params.pop_size} 0 0 \
          --release {params.pop_size} \
          --hom-antidote-effect {wildcards.hom} \
          --het-antidote-effect {wildcards.het} \
          --late-releases {params.late_count} \
          --late-releases-start {params.late_start} \
          --late-anti {params.late_anti} \
          --repetitions {params.repetitions} > {output}
      '''

rule all:
  input:
    rules.baseline.output,
    rules.simulations.input,
    rules.simulations_no_antidote_cost.input,
    rules.simulations_10.input,
    rules.simulations_10_no_antidote_cost.input,
    rules.simulations_late_wt.input,
    rules.simulations_late_wt_no_antidote_cost.input,
    rules.simulations_late_anti.input,
    rules.simulations_late_anti_no_antidote_cost.input,
    rules.simulations_late_anti_once.input,
    rules.simulations_late_anti_once_no_antidote_cost.input,
    rules.no_antidote.output.drive_25,
    rules.no_antidote.output.drive_50
