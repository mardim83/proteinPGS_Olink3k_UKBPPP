#!/bin/bash
#BSUB -q medium
#BSUB -n 5
#BSUB -J -fix_var_ids-
#BSUB -R "span[ptile=4]"
#BSUB -M 100000
module load ib
module load R/4.1.2
Rscript /hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/fix_var_ids_v3.R

