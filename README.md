large-cage
====

Agent-based modelling for the spread of gene drives and anti-drives in a large cage setting

Prerequisites
----

* python 3+
* numpy
* scipy
* scikit-learn
* snakemake

The necessary packages can be installed through `conda`:

    conda create -n large-cage numpy scipy scikit-learn snakemake
    conda activate large-cage

Or through `pip` using the system's python interpreter:

    python3 -m pip install numpy scipy scikit-learn snakemake

Usage
----

All analysis can be run using `snakemake`:

    snakemake -p all --cores CORES

where `CORES` is the number of cores used for parallelization.

A script to run a simulation in which drive and antidote heterozygous males
are introduced is provided (`src/simulation.py`), and can be run as follows:

    python3 src/simulation.py --parameters parameters/large/base/antidote.yaml

The `parameters/large/base/antidote.yaml` configuration file contains all the values
for the model parameters. The `parameters` folder contains all variations that were
used in the work.

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

The default parameters are defined on top of the `src/large_cage/agent.py` file, and can be overriden by providing a YAML file,
as shown in the example above.
