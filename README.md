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
