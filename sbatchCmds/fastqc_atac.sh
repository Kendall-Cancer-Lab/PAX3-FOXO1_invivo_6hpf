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
#SBATCH --array=0-7
#SBATCH --time=1-00:00:00

set -e ### stops bash script if line ends with error

echo ${HOSTNAME} ${SLURM_ARRAY_TASK_ID}

module purge
module load FastQC/0.11.9-Java-11.0.2

inputPath=/home/gdkendalllab/lab/raw_data/fastq/2023_04_03/JPK*

fileArray=(${inputPath}/*fastq.gz)

base_name=${fileArray[${SLURM_ARRAY_TASK_ID}]##*/}

fastqc \
  -o output/fastqc/atac \
  -t 2 \
  --extract \
  ${fileArray[${SLURM_ARRAY_TASK_ID}]}

rm output/fastqc/atac/${base_name%.fastq.gz}_fastqc.zip

perl scripts/splitFastqcData.pl \
    --input output/fastqc/atac/${base_name%.fastq.gz}_fastqc/fastqc_data.txt
