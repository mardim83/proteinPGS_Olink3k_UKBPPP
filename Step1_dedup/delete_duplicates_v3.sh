#!/bin/bash
#BSUB -q medium
#BSUB -n 5
#BSUB -J -de_dup-
#BSUB -R "span[ptile=4]"
#BSUB -M 100000
/lustre/workspace/projects/Genetics/bin/plink2 \
--pgen /lustre/workspace/projects/kimh132/resources/UKB/genotypes/UKB_500k_imputed/UKB_500k_imputed_Broad_EUR_x_Olink.pgen \
--psam /lustre/workspace/projects/kimh132/resources/UKB/genotypes/UKB_500k_imputed/UKB_500k_imputed_Broad_EUR_x_Olink.psam \
--pvar /hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/UKB_500k_imputed_Broad_EUR_x_Olink_marios.pvar \
--exclude /hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/UKB_500k_imputed_Broad_EUR_x_Olink_duplicates.txt \
--make-pgen \
--out /hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/UKB_500k_imputed_Broad_EUR_x_Olink_dedup

