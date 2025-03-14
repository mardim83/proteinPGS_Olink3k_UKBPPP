#!/bin/bash
#BSUB -q medium
#BSUB -n 5
#BSUB -J -identify_snps-
#BSUB -R "span[ptile=4]"
#BSUB -M 50000
module load ib
module load R/4.1.2
Rscript /hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/identify_snps_v3.R

