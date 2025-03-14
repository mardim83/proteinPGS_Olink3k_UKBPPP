#!/bin/bash
#BSUB -q medium
#BSUB -n 5
#BSUB -J -extract_QTLs_v3-
#BSUB -R "span[ptile=4]"
#BSUB -M 100000
#BSUB -o /hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/logs/extract_QTLs/extract_QTLs_v3.j.o
#BSUB -e /hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/logs/extract_QTLs/extract_QTLs_v3.j.e
module load ib
module load R/4.1.2
Rscript /hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/extract_QTLs_v3.R
