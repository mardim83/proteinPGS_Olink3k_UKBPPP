#!/bin/bash
#BSUB -q medium
#BSUB -n 5
#BSUB -J "extract_model_marios"
#BSUB -R "span[ptile=4]"
#BSUB -M 50000
#BSUB -o /hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/logs/genetic_score_training/v3/extract_model/extract_model_v3_marios.j.o
#BSUB -e /hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/logs/genetic_score_training/v3/extract_model/extract_model_v3_marios.j.e
module load ib
module load python/3.7.12
python3 /hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/extact_model_v3_marios.py 0 2601

