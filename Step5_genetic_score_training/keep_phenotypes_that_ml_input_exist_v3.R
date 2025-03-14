library(data.table)
olink_info= fread("/hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/ml_inputs/v3/phenotype_info.txt")
setwd("/hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/ml_inputs")
getwd()
myfolder = "v3"
allfiles = list.files(path=myfolder, pattern="*_variant_effects.txt", full.names=TRUE)
allfiles
allfiles = gsub("v3/|\\_variant_effects.txt$", "", allfiles)
ml_input_exist = data.table(allfiles)
colnames(ml_input_exist)[1] = "phenotype"
#olink_info_final = olink_info[phenotype %chin% unique(ml_input_exist$phenotype)]
#rm(olink_info_final)
olink_info_final=merge(ml_input_exist, olink_info, by="phenotype", all.x=T)
sum(duplicated(olink_info_final$phenotype))
library(tidyverse)
olink_info_final <- olink_info_final %>% distinct(phenotype, .keep_all = TRUE)
fwrite(olink_info_final[, .(phenotype,  UniProt, HGNC.symbol, chr, gene_start, Panel, UKBPPP_ProteinID)],
       sep="\t", quote=FALSE, file="/hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/ml_inputs/v3/phenotype_info_final.txt")
olink_info_final= fread("/hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/ml_inputs/v3/phenotype_info_final.txt")



library(data.table)
olink_pheno= fread("/hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/ml_inputs/v3/phenotypes.txt")
olink_pheno = olink_pheno[phenotype %chin% unique(olink_info_final$phenotype)]
fwrite(olink_pheno, sep="\t", quote=FALSE, file="/hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/ml_inputs/v3/phenotypes_final.txt")

#test = olink_pheno[phenotype %chin% "AAMDC_Q9H7C9_OID30236_v1_Cardiometabolic_II"]
