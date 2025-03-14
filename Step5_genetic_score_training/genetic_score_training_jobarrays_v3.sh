#!/bin/bash
#BSUB -q medium
#BSUB -n 4
#BSUB -J "genetic_score_training[0-2601]%10"
#BSUB -R "span[ptile=4]"
#BSUB -M 20000
#BSUB -o /hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/logs/genetic_score_training/v3/second_try/score_training_index_v3_%I.j.o
#BSUB -e /hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/logs/genetic_score_training/v3/second_try/score_training_index_v3_%I.j.e
module load ib
module load python/3.7.12
python3 /hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/run_olink3k_pgs_training_with_br_all_v3.py ${LSB_JOBINDEX} 0.000001 0.000001 0.000001 0.000001
