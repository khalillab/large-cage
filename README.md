large-cage
====

Agent-based modelling for the spread of gene drives and anti-drives in a large cage setting

Prerequisites
----

* python 3+
* numpy
* scipy
* snakemake

The necessary packages can installed through `conda`:

    conda create -n large-cage numpy scipy snakemake
    conda activate large-cage

Or through `pip` using the system's python interpreter:

    python3 -m pip install numpy scipy snakemake

Usage
----

All analysis can be run using `snakemake`:

    snakemake -p all --cores CPU

where `CPU` is the number of cores used for parallelization.

A script to run a simulation in which drive and antidote heterozygous males
is provided (`src/simulation.py`), and can be run as follows:

    python3 src/simulation.py --drive 71 72 --antidote 0 0 --pop-size 400 --repetitions 100 --time-points data/reference_time_points.txt

In the above example we are running 100 simulations with an initial releases
of 200 wild-type males, 200 wild-type females and 71 het. drive males. The second
introdcution is the same as the first one except that 72 het. drive individuals
are introduced on top of the 400 wild-type ones. The output is computed only
for those time points present in the `data/reference_time_points.txt` file
(one line per time point).

    python3 src/simulation.py --drive 71 72 100 --antidote 0 10 50 --pop-size 500 --repetitions 10 

In the second example we are running 10 simulations with three introductions:
on top of the 500 wild-type individuals we are introducing het. drive males
and het. antidote males, as follows:

* initial population: 71 drives, 0 antidotes
* first introduction: 72 drives, 10 antidotes
* second introduction: 100 drives, 50 antidotes

Since we didn't provide the `--time-points` argument, the script will report
output for each day of the simulations.

Output
----

The output is a tab-delimited table printed to `stdout` (which can be saved
to file through redirection, *e.g.* `> output.tsv`). The columns are the
following:

* round: identifier of the simulation
* time: time in days
* initial_release: wether at this time point one of the initial populations has been released
* pop: number of individuals inside the cage (adults + pupae)
* fpop: female individuals inside the cage (adults + pupae)
* eggs: number of eggs produced at this time point
* feggs: number of female eggs produced at this time point
* output: number of larvae + pupae that have developed outside the cage at this time point
* foutput: number of female larvae + pupae that have developed outside the cage at this time point
* fitness: proportion of female individuals that can mate at this time point
* WT: proportion of wild-type individuals inside the cage (*i.e.* genotype is WWWW)
* transgenes: proportion of non wild-type individuals inside the cage (*i.e.* any genotype that is not WWWW)
* drives: proportion of individuals with at least one drive allele inside the cage (*i.e.* any genotype containing D)
* antidote: proportion of individuals with at least one antidote allele inside the cage (*i.e.* any genotype containing A)
* resistance: proportion of individuals with at least one resistance allele inside the cage (*i.e.* any genotype containing R)
* DDAA, DDAW, DDWW, DRAA, DRAW, DRWW, DWAA, DWAW, DWWW, RRAA, RRAW, RRWW, RWAA, RWAW, RWWW, WWAA, WWAW, WWWW: proportion of individuals with each genotype
* DDAA.eggs, DDAW.eggs, DDWW.eggs, DRAA.eggs, DRAW.eggs, DRWW.eggs, DWAA.eggs, DWAW.eggs, DWWW.eggs, RRAA.eggs, RRAW.eggs, RRWW.eggs, RWAA.eggs, RWAW.eggs, RWWW.eggs, WWAA.eggs, WWAW.eggs, WWWW.eggs: same as above but for the eggs produced at this time point
* WT.output, transgenes.output, drives.output, antidote.output, resistance.output: same as the previous columns, but for the larvae + pupae maturing outside the cage
* DDAA.output, DDAW.output, DDWW.output, DRAA.output, DRAW.output, DRWW.output, DWAA.output, DWAW.output, DWWW.output, RRAA.output, RRAW.output, RRWW.output, RWAA.output, RWAW.output, RWWW.output, WWAA.output, WWAW.output, WWWW.output: proportion of each individual genotype for the larvae + pupae maturing outside the cage
