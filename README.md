# sourmash-oddify

sourmash-oddify contains tools and workflows for examining bacterial
and archaeal genome contamination with k-mers using
[sourmash](http://sourmash.rtfd.io/).

Briefly, sourmash-oddify will:

* classify genomes with GTDB-Tk
* build sourmash LCA databases based on the resulting taxonomy
* identify, align, report, and ~remove taxonomically discordant sub-sequences

The primary value-add of sourmash-oddify is that it uses sourmash to
quickly identify genomes that share many taxonomically discordant
k-mers. These genomes can then be further investigated.  (See the
Kraken section in [this blog
post](http://ivory.idyll.org/blog/2017-something-about-kmers.html) for
more info on how we identify taxonomically discordant k-mers.)

sourmash-oddify is in the early stages of development but we are happy
to support your using it! Please post questions and problems to [our
issue tracker](https://github.com/dib-lab/sourmash-oddify/issues).

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
More on what this output is ...later. :)

### Configuring prefixes

For now, your genome identifiers must have a common prefix of at least one
character. Put this prefix in the `filter_prefixes` in your config file.

### GTDB-Tk reference data

Since sourmash-oddify runs gtdbtk, you'll need the [GTDB-Tk reference
data](https://github.com/Ecogenomics/GTDBTk/blob/stable/README.md).
Download it, unpack it somewhere, and then set the value for
`gtdbtk_data` in your config file to the unpacked directory path.
Note, this should be set to the top level directory containing the
subdirectories `fastani/` and `taxonomy/`, e.g.  `release89/`.

## Acknowledgements

sourmash-oddify relies on the [GTDB
taxonomy](https://www.biorxiv.org/content/10.1101/256800v2) and the
[GTDB-Tk toolkit](https://github.com/Ecogenomics/GtdbTk). Thanks, nice
people!

