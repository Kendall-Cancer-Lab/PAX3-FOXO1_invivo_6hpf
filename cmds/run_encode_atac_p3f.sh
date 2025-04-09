#!/bin/bash

#SBATCH --job-name=run_encode_atac

#SBATCH --time=12:00:00

#SBATCH --partition=general

#SBATCH --output=run_encode_atac_p3f_%j.out

#SBATCH --account=gdkendalllab

#ml ENCODE/caper/2.1.0
ml ENCODE/caper/2.3.2
ml Java/18 
set -x
INPUT_JSON="pax3foxo1_atac.json"
caper hpc submit atac.wdl -i "${INPUT_JSON}" --singularity \
	--leader-job-name p3f_atac_encode
