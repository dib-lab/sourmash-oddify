#!/bin/bash -login
#SBATCH -p bmm
#SBATCH -J run-gtdbtk
#SBATCH --mail-type=ALL
#SBATCH --mail-user=titus@idyll.org
#SBATCH -t 3-0:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=200gb

. "/home/ctbrown/miniconda3/etc/profile.d/conda.sh"

cd /home/ctbrown/2019-sourmash-gtdb/gtdbtk

conda activate sgc

snakemake --use-conda -j 8

set -o nounset
set -o errexit
set -x

env | grep SLURM            # Print out values of the current jobs SLURM environment variables

scontrol show job ${SLURM_JOB_ID}     # Print out final statistics about resource uses before job exits

