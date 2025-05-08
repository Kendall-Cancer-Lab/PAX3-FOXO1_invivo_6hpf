#!/bin/bash

#SBATCH --job-name=r_homer_motif

#SBATCH --time=12:00:00

#SBATCH --partition=general

#SBATCH --output=r_homer_motif_%j.out

#SBATCH --account=gdkendalllab

#SBATCH --cpus-per-task=10

ml homer/4.11.1
set -x
echo $peakfile
echo $outdir
findMotifsGenome.pl \
	$peakfile \
	/gpfs0/home/gdkendalllab/lab/references/fasta/danRer11.fa \
	$outdir \
	-size given \
	-p $SLURM_CPUS_PER_TASK
