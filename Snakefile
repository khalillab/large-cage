rule files:
  params:
    time_points = 'data/reference_time_points.txt'

files = rules.files.params

rule no_antidote:
  message:
    '''
    Run simulations without the antidote
    '''
  input:
      time_points = files.time_points
  output:
      drive_25 = 'out/no_antidote/drive_25.tsv',
      drive_50 = 'out/no_antidote/drive_50.tsv'
  params:
      pop_size = 400,
      repetitions = 100,
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
          --drive-1 {params.drive_50_1} {params.drive_50_2} \
          --antidote 0 0 \
          --pop-size {params.pop_size} \
          --repetitions {params.repetitions} \
          --time-points {input.time_points} > {output.drive_50}
      '''

rule all:
  input:
    rules.no_antidote.output.drive_25,
    rules.no_antidote.output.drive_50
