# --------------------------------------------------------------------------------------
# Load libraries/dependencies
# --------------------------------------------------------------------------------------
library(tidyverse)
library(data.table)
library(openxlsx)
library(foreach)
library(doMC)

# --------------------------------------------------------------------------------------
# Load in the ldthinned data
# ---------------------------------------------------------------------------------------
varset = fread("/hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/UKB_500k_imputed_Broad_EUR_x_Olink_dedup_unambig_SNPs_maf0.005_ldthin0.8.pvar")
setnames(varset, c("chr", "pos", "id", "ref", "alt"))
varset$pos = as.character(varset$pos)
varset$chr_pos = paste(varset$chr, varset$pos, sep=":")

# --------------------------------------------------------------------------------------
# Load in the olink data
# ---------------------------------------------------------------------------------------
#-log10(0.00000005)  # [1] 7.30103 that will be used as log10pval cuttoff in later stage
#-log10(0.0001)  # [1] 4 that will be used as log10pval cuttoff to ensure files will be writen for all proteins

proteins = "/hpc/grid/hgcb/workspace/projects/P181_pQTLs_UKB_Olink/summary_stats_discovery/v2"
protein_files = list.files(path=proteins)

protein_files <- protein_files[! protein_files %in% c('UKBPPP_GWAS_v2_RS2_manifest.txt')]   
#protein_files <- c("ABHD14B_Q96IU4_OID20921_v1_Neurology", "ABL1_P00519_OID21280_v1_Oncology")
#olink_ss = fread("/hpc/grid/hgcb/workspace/projects/P181_pQTLs_UKB_Olink/summary_stats_discovery/v2/ABHD14B_Q96IU4_OID20921_v1_Neurology/discovery_RS2_chunk2_chr1_ABHD14B:Q96IU4:OID20921:v1:Neurology.regenie.gz")

#protein_files <- head(protein_files, n=1500)
#protein_files <- tail(protein_files, n=1440)
#protein_files = c("A1BG_P04217_OID30771_v1_Inflammation_II")

for(i in {protein_files})  {
olink_GWAS = paste("/hpc/grid/hgcb/workspace/projects/P181_pQTLs_UKB_Olink/summary_stats_discovery/v2", i, sep='/')
  
olink_files = list.files(path=olink_GWAS, pattern="*.gz$")
olink_ss = foreach(ff = olink_files, .combine=rbind) %do% {
  ss = fread(cmd = sprintf('zcat %s/%s | tail -n +2 | awk \'BEGIN { OFS="\t" } { if ( $13 > 4 ) { print $3, $5, $4, $10, $13 } }\'', olink_GWAS, ff))
  setnames(ss, c("ID","effect_allele", "other_allele", "beta", "log10pval"))
  ss$chr_pos = stringr::str_extract(ss$ID, "[^:]*:[^:]*")
  ss = ss[varset[, .(chr_pos)], on = .(chr_pos), nomatch=0]
  return(ss)
  }
olink_ss = olink_ss %>% separate(ID, c('chr', 'pos', "A0", "A1", "method", "version"), sep = ":")
olink_ss = subset(olink_ss, select = -c(A0, A1, method, version, chr_pos))
olink_ss$phenotype = gsub(".regenie.gz", "", gsub("_chr*", "", ff))
fwrite(olink_ss, file=paste("/hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/ldthinned/collate_sumstats_discovery_v3/", i, "_olink_log10p_4.txt",sep=''), sep="\t", quote=F)
rm(olink_ss)
gc()
}

setwd("/hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/ldthinned/")
getwd()

myfolder = "collate_sumstats_discovery_v3"
allfiles = list.files(path=myfolder, pattern="*.txt", full.names=TRUE)
allfiles

data <- lapply(allfiles, read.table, sep="\t", header=TRUE)
data <- do.call(rbind, data)

fwrite(data, file="/hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/ldthinned/sumstats_discovery/olink_log10p_4_v3.txt", sep="\t", quote=F)

