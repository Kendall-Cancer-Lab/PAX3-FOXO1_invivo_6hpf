#!/bin/bash

#SBATCH --job-name=run_samtools_merge

#SBATCH --partition=general

module load SAMtools/1.15
set -x
samtools merge - $infiles | samtools sort -o $outfile -
samtools index $outfile
