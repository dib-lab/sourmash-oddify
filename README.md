# sourmash-oddify

Tools and workflows for examining bacterial and archaeal genome contamination
with k-mers a la [sourmash](http://sourmash.rtfd.io/.

## Quickstart

You can install all of the software dependencies with
[bioconda](https://bioconda.github.io/) by running:

```
conda env create -f environment.yml
```

The primary way to use sourmash-oddify is to run the snakemake
workflow in the Snakefile. Do so by:

1. copying `conf/default.yml` to a new file and editing it
2. running `snakemake --use-conda --configfile=newconf.yml`

The output will be placed in the `outputs_dir` specified in your config file.
More on what it is ...later. :)

## Other dependencies

sourmash-oddify relies on the [GTDB
taxonomy](https://www.biorxiv.org/content/10.1101/256800v2) and the
[GTDB-Tk toolkit](https://github.com/Ecogenomics/GtdbTk).

TODO: add instructions re @@download gtdbtk database; release89

### Setting GTDBTK_DATA_PATH

One tricky bit is the need to set the GTDBTK_DATA_PATH variable
properly, so that the gtdbtk pipeline can find the GTDB release.  For
now, this needs to be done manually in a file that has to be created
by snakemake.

First, run `snakemake --use-conda` to set up the gtdbtk environment.
(It will end in an error because `GTDBTK_DATA_PATH` is not set!)

Then, find the environment activation script named `gtdbtk.sh`; you can
use the following command to find the file.

`ls $(find .snakemake/conda -name activate.d)/gtdbtk.sh`

Edit the file so that `GTDBTK_DATA_PATH` is set to the top level directory
in release89 containing the subdirectories `fastani/` and `taxonomy/`.
