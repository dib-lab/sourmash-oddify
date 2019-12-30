# sourmash-oddify

Tools and workflows for examining bacterial and archaeal genome contamination
with k-mers a la [sourmash](http://sourmash.rtfd.io/.

## Basic technical overview.

You'll need [snakemake](https://snakemake.readthedocs.io/en/stable/)
and [bioconda](https://bioconda.github.io/).

The primary way to use sourmash-oddify is to run the snakemake
workflow in the Snakefile, and to configure it by

## Other dependencies

sourmash-oddify relies on the [GTDB
taxonomy](https://www.biorxiv.org/content/10.1101/256800v2) and the
[GTDB-Tk toolkit](https://github.com/Ecogenomics/GtdbTk).

## Installation

@@download gtdbtk database; release89

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
