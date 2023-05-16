#!/bin/sh
#SBATCH --account=gdkendalllab
#SBATCH --error=slurmOut/fastqc-%j.txt
#SBATCH --output=slurmOut/fastqc-%j.txt
#SBATCH --mem=2G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --job-name fastqc
#SBATCH --wait
#SBATCH --array=0-63
#SBATCH --time=1-00:00:00

set -e ### stops bash script if line ends with error

echo ${HOSTNAME} ${SLURM_ARRAY_TASK_ID}

module purge
module load FastQC/0.11.9-Java-11.0.2

inputPath=/home/gdkendalllab/lab/raw_data/fastq/2023_03_17/JPK*

fileArray=(${inputPath}/*fastq.gz)

base_name=${fileArray[${SLURM_ARRAY_TASK_ID}]##*/}

fastqc \
  -o output/fastqc \
  -t 2 \
  --extract \
  ${fileArray[${SLURM_ARRAY_TASK_ID}]}

rm output/fastqc/${base_name%.fastq.gz}_fastqc.zip

perl scripts/splitFastqcData.pl \
    --input output/fastqc/${base_name%.fastq.gz}_fastqc/fastqc_data.txt \
    --outputDir output/fastqc/${base_name%.fastq.gz}_fastqc/