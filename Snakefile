configfile: 'config.yml'

filter_prefixes = config['filter_prefixes']
genomes_extension = config['genomes_extension']
scaled = int(config['scaled'])
build_ksizes = [ int(k) for k in config['build_ksizes'] ]
outputs_dir = config.get('outputs_dir', 'outputs')
genomes_dir = config.get('genomes_dir', 'genomes')

###

require_taxonomy_bool=0

###

import os
import random

all_files = []
for root, dirs, files in os.walk(genomes_dir, topdown=False):
    for name in files:
        if name.endswith(genomes_extension):
            filename = os.path.join(root, name)
            if filename.startswith('./'): filename = filename[2:]
            all_files.append(filename)

rule all:
    input:
        os.path.join(outputs_dir, "gtdbtk/"),
        [ i + '.sig' for i in all_files ],
        expand("{outprefix}/{prefix}-genomes-k{k}.lca.json.gz", k=build_ksizes,
               prefix=filter_prefixes, outprefix=outputs_dir),
        expand("{outprefix}/{prefix}-oddities-k{k}.txt", k=build_ksizes,
                prefix=filter_prefixes, outprefix=outputs_dir),
        expand("{outprefix}/{prefix}-oddities-k{k}.examine.txt",
               k=build_ksizes, prefix=filter_prefixes, outprefix=outputs_dir)

rule sigs:
    input:
        "{genome}"
    output:
        "{genome}.sig"
    params:
        extension=genomes_extension
    shell:
        "sourmash compute -k 21,31,51 --scaled=1000 {input} -o {output} --merge=$(basename {input} {params.extension})"

rule gtdbtk_gather_matches:
    """
    this rule require the gtdbtk databases. The tool finds the database by 
    using a path specified in a file in the environment. I predownloaded the 
    databases and placed them in the required location.
    The path is in this file:
    .snakemake/conda/1261315d/etc/conda/activate.d/gtdbtk.sh
    """
    input: genomes_dir
    output: directory(expand("{outprefix}/gtdbtk/", outprefix=outputs_dir))
    params:
        outdir = expand("{outprefix}/gtdbtk", outprefix=outputs_dir),
        extension=genomes_extension
    conda: "env-gtdbtk.yml"
    shell:'''
    gtdbtk classify_wf --genome_dir {input} --out_dir {params.outdir} --cpus 8 --extension {params.extension}
    '''

rule make_lineages_csv:
    input:
        expand("{outprefix}/gtdbtk/", outprefix=outputs_dir)
    output:
        expand("{outprefix}/{{filter_prefix}}-lineages.csv", outprefix=outputs_dir)
    params:
        outputs_dir=outputs_dir
    shell:
        "../gtdbtk-to-lineages-csv.py {params.outputs_dir}/gtdbtk/ {output} --filter-prefix={wildcards.filter_prefix:q}"

rule make_lca_db:
    input:
        os.path.join(outputs_dir, "{prefix}-lineages.csv"),
    output:
        os.path.join(outputs_dir, "{prefix}-genomes-k{ksize}.lca.json.gz")
    params:
        scaled=scaled,
        require_taxonomy_arg="--require-taxonomy",
        genomes_dir=genomes_dir
    shell:
        "sourmash lca index {params.require_taxonomy_arg} {input} {output} --traverse-directory {params.genomes_dir} -k {wildcards.ksize} --scaled={params.scaled}"

rule make_oddities_txt:
    input:
        os.path.join(outputs_dir, "{prefix}-genomes-k{ksize}.lca.json.gz")
    output:
        os.path.join(outputs_dir, "{prefix}-oddities-k{ksize,[0-9]+}.txt"),
        os.path.join(outputs_dir, "{prefix}-oddities-k{ksize,[0-9]+}.csv")
    params:
        outputs_dir=outputs_dir
    shell:
        "../find-oddities.py {input} --lowest-rank=superkingdom --minimum-hashes=10 --prefix={params.outputs_dir}/{wildcards.prefix}-oddities-k{wildcards.ksize} > {output[0]}"

rule make_oddities_examine_txt:
    input:
        os.path.join(outputs_dir, "{prefix}-oddities-k{ksize}.txt"),
        os.path.join(outputs_dir, "{prefix}-genomes-k{ksize}.lca.json.gz"),
        os.path.join(outputs_dir, "{prefix}-oddities-k{ksize}.csv")
    output:
        os.path.join(outputs_dir, "{prefix}-oddities-k{ksize,[0-9]+}.examine.txt")
    params:
        outputs_dir=outputs_dir,
        genomes_dir=genomes_dir
    shell:
        "../find-oddities-examine.py {params.outputs_dir}/{wildcards.prefix}-oddities-k{wildcards.ksize}.csv {params.genomes_dir} --percent-threshold=95 --length-threshold=0 > {output}"
