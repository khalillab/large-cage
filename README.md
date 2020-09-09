large-cage
====

Agent-based modelling for the spread of gene drives and anti-drives in a large cage setting

Prerequisites
----

* python 3+
* numpy
* scipy
* snakemake

The necessary packages can be installed through `conda`:

    conda create -n large-cage numpy scipy snakemake
    conda activate large-cage

A version freeze of the conda environment is also provided in the `conda.yaml`
file and can be used to install the requirements:

    conda env create --file conda.yaml
    conda activate large-cage

Or through `pip` using the system's python interpreter:

    python3 -m pip install numpy scipy snakemake

Usage
----

All analysis can be run using `snakemake`:

    snakemake -p all --cores CPU

where `CPU` is the number of cores used for parallelization.

A script to run a simulation in which drive and antidote heterozygous males
are introduced is provided (`src/simulation.py`), and can be run as follows:

    python3 src/simulation.py --wild-type 400 400 --drive 71 72 --antidote 0 0 --release 400 --repetitions 100 --time-points data/reference_time_points.txt

In the above example we are running 100 simulations with an initial release
of 200 wild-type males, 200 wild-type females and 71 het. drive males. The second
introduction is the same as the first one except that 72 het. drive individuals
are introduced on top of the 400 wild-type ones. The output is computed only
for those time points present in the `data/reference_time_points.txt` file
(one line per time point).

    python3 src/simulation.py --wild-type 400 400 400 --drive 71 72 100 --antidote 0 10 50 --release 500 --repetitions 10 

In the second example we are running 10 simulations with three introductions:
on top of the 400 wild-type individuals we are introducing het. drive males
and het. antidote males, as follows:

* initial population: 71 drives, 0 antidotes
* first introduction: 72 drives, 10 antidotes
* second introduction: 100 drives, 50 antidotes

The `--release` option indicates how many pupae are introduced once the eggs produced in the 
cage are developed.
Since we didn't provide the `--time-points` argument, the script will report
output for each day of the simulations.

Additional releases
----

The `simulation.py` script has the option to perform additional introductions
of either wild-type males or heterozygous antidote males. The three following
arguments have to be provided:

* `--late-releases`: number of extra releases (-1 indicates continuous releases)
* `--late-releases-start`: Start time to start the extra releases
* `--late-antidote` and `--late-wt`: how many males to release (antidote or WT, respectively)

Example:

    python3 src/simulation.py --drive 100 100 --antidote 400 400 --wild-type 400 400 --late-releases -1 --late-releases-start 30 --late-wt 100

Will start releasing 100 male WT individuals starting from day 30 at each feeding cycle, until the end of the simulation.

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

Changing parameters
----

The default parameters are defined on top of the `src/large_cage/agent.py` file, and can be overriden in two ways:

* From the script launching the simulation (see the `src/simulation.py` for an example, only a few parameters)
* By changing the `parameters.py` file and using the `--override-parameters` option of the `src/simulation.py` script

It is NOT advisable to change the parameters inside the
`src/large_cage/agent.py` file manually, unless the changes
are being committed to version control.
Also note that if the `--override-parameters` option is used, then the parameters written in the `src/parameters.py` file have priority over
those provided in the command line.
