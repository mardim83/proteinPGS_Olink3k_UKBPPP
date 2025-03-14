library(data.table)

# Remove extra row identifier tacked onto the variant IDs
  pvar = fread("/hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/UKB_500k_imputed_Broad_EUR_x_Olink_dedup.pvar")
  pvar[, ID := gsub(":[0-9]+?$", "", ID)]
  fwrite(pvar, sep="\t", quote=FALSE, file="/hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/UKB_500k_imputed_Broad_EUR_x_Olink_dedup.pvar")