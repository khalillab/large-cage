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
      drive_25_1 = 71,
      drive_25_2 = 72,
      drive_50_1 = 142,
      drive_50_2 = 143
  shell:
      '''
      python3 src/simulation_noantidote.py \
          --drive-1 {params.drive_25_1} \
          --drive-2 {params.drive_25_2} \
          --time-points {input.time_points} > {output.drive_25}
      python3 src/simulation_noantidote.py \
          --drive-1 {params.drive_50_1} \
          --drive-2 {params.drive_50_2} \
          --time-points {input.time_points} > {output.drive_50}
      '''

rule all:
  input:
    rules.no_antidote.output.drive_25,
    rules.no_antidote.output.drive_50
