#!/bin/sh
#SBATCH --account=gdkendalllab
#SBATCH --array=0-7
#SBATCH --error=slurmOut/macs2_callPeaks2-%j.txt
#SBATCH --output=slurmOut/macs2_callPeaks2-%j.txt
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --job-name macs2Peakcall
#SBATCH --wait
#SBATCH --time=2-08:00:00

set -e ### stops bash script if line ends with error

echo $HOSTNAME $SLURM_ARRAY_TASK_ID

module purge
module load GCC/7.3.0-2.30 \
            OpenMPI/3.1.1 \
            MACS2/2.2.5-Python-3.6.6

fileArray=(JPK0028
           JPK0030
           JPK0032
           JPK0034
           JPK0036
           JPK0038
           JPK0040
           JPK0042)

controlArray=(JPK0029
              JPK0031
              JPK0033
              JPK0035
              JPK0037
              JPK0039
              JPK0041
              JPK0043)

echo ${fileArray[${SLURM_ARRAY_TASK_ID}]

input_path=/home/gdkendalllab/lab/analyses/rsjxk002/P3F_ChIP_6hpf/align_nosplice

macs2 callpeak \
  -t output/align/chip/${fileArray[${SLURM_ARRAY_TASK_ID}]}_dedupd.bam \
  -c output/align/chip/${controlArray[${SLURM_ARRAY_TASK_ID}]}_dedupd.bam \
  --outdir output/peaks \
  -n ${fileArray[${SLURM_ARRAY_TASK_ID}]} \
  -f BAMPE \
  -B \
  -g 1679203469 \
  --keep-dup all
