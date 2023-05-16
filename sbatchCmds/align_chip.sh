#!/bin/sh
#SBATCH --error=slurmOut/align-%j.txt
#SBATCH --output=slurmOut/align-%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --job-name alignChIP
#SBATCH --partition=himem,general
#SBATCH --array=0-15
#SBATCH --wait
#SBATCH --time=1-00:00:00

set -e

echo ${HOSTNAME} ${SLURM_ARRAY_TASK_ID}

module purge
module load GCC/9.3.0 \
            GCCcore/9.3.0 \
            HISAT2/2.2.1 \
            SAMtools/1.15

fileArray=($(ls /home/gdkendalllab/lab/raw_data/fastq/2023_03_17/JPK*/*_L001_R1_001.fastq.gz))

R1=${fileArray[$SLURM_ARRAY_TASK_ID]}
R2=${R1/R1_001.fastq.gz/R2_001.fastq.gz}

echo $R1

baseName=${R1##*/}
baseName=${baseName%%_*}

hisat2 \
        -x /home/gdkendalllab/lab/references/hisat2/danRer11 \
        -1 ${R1},${R1/_L001_/_L002_} \
        -2 ${R2},${R2/_L001_/_L002_} \
        --threads 15 \
        -k 1 \
        --no-spliced-alignment \
        --no-temp-splicesite \
        --no-mixed \
        --no-discordant \
        --summary-file ${baseName}Summary.txt \
    | samtools fixmate \
        -@ 10 \
        -m \
        - - \
    | samtools sort \
        -T /gpfs0/scratch/junkdir1337/ \
        -@ 10 \
        -m 15G \
        -O BAM \
    | samtools view \
        -@ 10 \
        -h \
        - \
    | python scripts/count_inline.py \
        -i "^@" \
        -o output/counts/chip_post_align_${baseName}.tsv \
    | samtools markdup \
        -@ 10 \
        -O SAM \
        -r \
        - - \
    | python scripts/count_inline.py \
        -i "^@" \
        -o output/counts/chip_post_rmdup_${baseName}.tsv \
    | samtools view \
        -@ 10 \
        -b \
        - \
    > output/align/chip/${baseName}_dedupd.bam

samtools index output/align/chip/${baseName}_dedupd.bam
