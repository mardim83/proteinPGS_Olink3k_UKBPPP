#!/bin/bash
#BSUB -q medium
#BSUB -n 5
#BSUB -J -collate_QTLs_v3-
#BSUB -R "span[ptile=4]"
#BSUB -M 100000
#BSUB -o /hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/sumstats_discovery/olink_log10p_4_v3.txt.j.o
#BSUB -e /hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/sumstats_discovery/olink_log10p_4_v3.txt.j.e
module load ib
module load R/4.1.2
Rscript /hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/collate_QTLs_v3.R
