#!/bin/bash
#BSUB -q medium
#BSUB -n 5
#BSUB -J -maf0.005_ldthin0.8_v3-
#BSUB -R "span[ptile=4]"
#BSUB -M 50000
/lustre/workspace/projects/Genetics/bin/plink2 \
--pfile /hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/UKB_500k_imputed_Broad_EUR_x_Olink_dedup \
--extract /hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/ldthinned/keep_v3.txt \
--maf 0.005 \
--indep-pairwise 1000kb 0.8 \
--silent \
--out /hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/ldthinned/ldthinned_v3

