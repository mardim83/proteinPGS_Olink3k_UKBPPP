# --------------------------------------------------------------------------------------
# Load libraries/dependencies
# --------------------------------------------------------------------------------------
library(data.table)
library(foreach)
library(doMC)
library(AnnotationHub) # Bioconductor package
library(annotables) # remotes::install_github("stephenturner/annotables")

#usethis::create_github_token()
#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
#usethis::edit_r_environ()


ncores = 25
ncores = 5

parallelise_fread = function() {
  setDTthreads(ncores)
  registerDoMC(1)
}

parallelise_foreach = function() {
  setDTthreads(1)
  registerDoMC(ncores)
}


## ======================================================================================
## First, we want to load the summary statistics for all platforms, and filter to the
## ld-thinned variant set, and add a basic filter of P < 0.01 - we want to find a
## reasonable P-value threshold across all platforms but need to load the summary stats
## for all measurements
## ======================================================================================


varset = fread("/hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/UKB_500k_imputed_Broad_EUR_x_Olink_dedup_unambig_SNPs_maf0.005_ldthin0.8.pvar")
setnames(varset, c("chr", "pos", "id", "ref", "alt"))
setkey(varset, chr, pos)


#work with one example here (THEN DELETE IN FINAL R script)
#olink_SS = fread("/hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/ldthinned/collate_sumstats_discovery/A1BG_P04217_OID30771_v1_Inflammation_II_olink_log10p_4.txt")
#ph1 = unique(olink_SS$phenotype)

#test = sub('.*?[_]', '', ph1)
#test = sub('.*?[_]', '', test)

#olink_SS$phenotype = sub('.*?[_]', '', olink_SS$phenotype)
#olink_SS$phenotype = sub('.*?[_]', '', olink_SS$phenotype)
#head(olink_SS$phenotype, n=200)


#THE ACTUAL MERGED FILE TO BE USED AT THE END IS THIS ONE:
"/hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/ldthinned/sumstats_discovery/olink_log10p_4_v3.txt"
olink_SS = fread("/hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/ldthinned/sumstats_discovery/olink_log10p_4_v3.txt")

# Make variable names match across data.tables
olink_SS$phenotype = sub('.*?[_]', '', olink_SS$phenotype)
olink_SS$phenotype = sub('.*?[_]', '', olink_SS$phenotype)


olink_pheno = fread("/hpc/grid/hgcb/workspace/projects/P015_UKBiobank_Resource/HLA_data_500K/input_files/discovery_pheno_forconsortium_v1.tsv")

olink_info = fread("/hpc/grid/wip_drm_targetsciences/projects/p070_CVITEN_Jasmine/UKB_OLINK_shared_core_files/V2_Expansion3072/olink_protein_map_3k_v1.tsv")
olink_info$phenotype = paste(olink_info$UKBPPP_ProteinID, olink_info$Panel, sep=":")
#successfull match of variable names between "olink_SS", "olink_pheno" and "olink_info"


key = fread("/hpc/grid/hgcb/workspace/projects/P015_UKBiobank_Resource/HLA_data_500K/input_files/combined_sample_map_v2.txt")
#olink_pheno_WRONG_OLD = merge(olink_pheno, key, by.x="pseudo_ind_id", by.y="pseudo_ind_id_OLD")
#olink_pheno_CORRECT_NEW = merge(olink_pheno, key, by.x="pseudo_ind_id", by.y="pseudo_ind_id_NEW")


#keyfile3 = fread("/hpc/grid/hgcb/workspace/projects/P015_UKBiobank_Resource/linking_helper_files/three.application.linking.file.allsubjects.txt")
#keyfile3$f.eid.AMRA = as.numeric(keyfile3$f.eid.AMRA)
#keyfile3$f.eid.Regn = as.numeric(keyfile3$f.eid.Regn)
#keyfile3$f.eid.Pfe = as.numeric(keyfile3$f.eid.Pfe)
#checkIID = fread("/lustre/workspace/projects/kimh132/resources/UKB/genotypes/UKB_500k_imputed/UKB_500k_imputed_Broad_EUR_x_Olink.psam")
#test1 = merge(checkIID, keyfile3, by.x="IID", by.y="f.eid.AMRA")
#test2 = merge(checkIID, keyfile3, by.x="IID", by.y="f.eid.Regn")
#test3 = merge(checkIID, keyfile3, by.x="IID", by.y="f.eid.Pfe")
#Perfect match between this new genotype file based IID with "f.eid.Regn" THEREFORE I NEED TO TRANSLATE THE pseudo_ind_id to f.eid.Regn and name as IID

library(dplyr) 
key = key %>% select(pseudo_ind_id_NEW, f.eid.Regn)
olink_pheno = merge(olink_pheno, key, by.x="pseudo_ind_id", by.y="pseudo_ind_id_NEW")
 
olink_pheno = olink_pheno %>% relocate(f.eid.Regn)
#I WILL RENAME "f.eid.Regn" to "IID"
colnames(olink_pheno)[1] = "IID"
olink_pheno$pseudo_ind_id = NULL
nrow(olink_pheno) # [1] 34557
sum(duplicated(olink_pheno))

checkIID = fread("/lustre/workspace/projects/kimh132/resources/UKB/genotypes/UKB_500k_imputed/UKB_500k_imputed_Broad_EUR_x_Olink.psam")
nrow(checkIID) # [1] 46571

test = merge(olink_pheno, checkIID, by = "IID")
nrow(test) # [1] 34490



# Now test load summary stats for all SNPs passing pthresh genome wide 
#Cuttoff to be applied
#-log10(0.00000001)  # [1] 8
#-log10(0.00000005)  # [1] 7.30103

#Initially test only for 1 random unique IDs "A1BG:P04217:OID30771:v1:Inflammation_II"
#this_ss = fread("/hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/ldthinned/collate_sumstats_discovery_v2/A1BG_P04217_OID30771_v1_Inflammation_II_olink_log10p_4.txt")
#this_ss = this_ss[, .(chr=chr, pos=pos, effect_allele=effect_allele, other_allele=other_allele, effect=beta, log10pval=log10pval)] 
#this_ss = this_ss[log10pval > 7.30103]

#olink_GWAS = "/hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/ldthinned/collate_sumstats_discovery_v2"
#var_id= unique(olink_SS$phenotype)
#var_id = gsub(":","_", olink_SS[phenotype == "A1BG:P04217:OID30771:v1:Inflammation_II", unique(phenotype)])



# --------------------------------------------------------------------------------------
# Set paths
# --------------------------------------------------------------------------------------

out_dir = "/hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/ml_inputs/v3"

tmpdir = sprintf("%s/tmpdir", out_dir)

olink_GWAS = "/hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/ldthinned/collate_sumstats_discovery_v3"
#test2IDs = c("ABHD14B:Q96IU4:OID20921:v1:Neurology", "ABL1:P00519:OID21280:v1:Oncology")
#HOWEVER IF THIS SUCCEED THEN CHANGE THE test2IDs with unique(olink_SS$phenotype)

parallelise_foreach()
foreach(phen_id = unique(olink_SS$phenotype), .combine=c) %dopar% {
  phen_ss = foreach(var_id = gsub(":","_", olink_SS[phenotype == phen_id, unique(phenotype)]), .combine=rbind) %do% {
    # Load full summary stats

    this_ss = fread(sprintf("%s/%s_olink_log10p_4.txt", olink_GWAS, var_id), tmpdir=tmpdir)
    this_ss$chr = as.character(this_ss$chr)
    this_ss = this_ss[, .(chr=chr, pos=pos, effect_allele=effect_allele, other_allele=other_allele, effect=beta, log10pval=log10pval)] 
    this_ss = this_ss[log10pval > 7.30103]
    
    # Filter to varset snps
    if (nrow(this_ss) > 0) {
      this_ss = rbind(
        this_ss[varset[, .(chr, pos, ref, alt)], on = .(chr, pos, effect_allele=alt, other_allele=ref), nomatch=0],
        this_ss[varset[, .(chr, pos, ref, alt)], on = .(chr, pos, effect_allele=ref, other_allele=alt), nomatch=0]
      )
      this_ss = this_ss[!is.na(effect)] 
    }
    return(this_ss)
  }
  
  
  # R, do you even garbage collect? WTF.
  if(exists("this_ss")) { rm(this_ss) }
  gc()
  
  # Remove variants that did not pass the p-value threshold for all panels for each protein
  npan = olink_info[phenotype == phen_id, length(unique(phenotype))]
  varn = phen_ss[,.N, by=.(chr, pos, effect_allele, other_allele)]
  pass = varn[N == npan]
  pass[,N := NULL]
  phen_ss = phen_ss[pass, on = .(chr, pos, effect_allele, other_allele)]

  
  # Get SNPs at genome-wide P < trans_pthresh 
  phen_ss = phen_ss[log10pval > 7.30103]
  
  # If any QTLs, proceed
  if (nrow(phen_ss) > 0) {
    
    phen_ss = phen_ss[, .SD[which.max(log10pval)], by=.(chr, pos)] # some duplicate results (duplicate SNPs with different INFO). Take best estimate (highest log10pval).
    phen_ss[varset, on = .(chr, pos), rsid := id] # add rsid
    phen_ss = phen_ss[, .(rsid, chr, pos, effect_allele, other_allele, effect, log10pval)][order(pos)][order(chr)]
    
    # write out
    fwrite(phen_ss, sep="\t", quote=FALSE, file=sprintf("%s/%s_variant_effects.txt", out_dir, var_id))
    
    # Free up memory and garbage collect
    rm(phen_ss)
    gc()
    
    return(phen_id)
  }
}
##DONE SUCCESSFULL



#Code for reshaping based pseudo_ind_id and rename to IID: 
olink_pheno = reshape2::melt(olink_pheno, id.vars = c("IID") )


# Filter phenotype data 
library(tidyverse)
olink_pheno = olink_pheno %>% drop_na(value)

#Filter info sheet and phenotype data to measurements with at least 1 variant passing the P-value threshold
str(olink_pheno)
olink_pheno$variable = as.character(olink_pheno$variable)
colnames(olink_pheno)[2] = "phenotype"
olink_pheno = as.data.table(olink_pheno)

olink_info = olink_info[phenotype %chin% unique(olink_SS$phenotype)]
#change again from : to _ in order to be easier to read in the next script 
olink_info$phenotype = gsub(":","_", olink_info$phenotype)
# write out 
fwrite(olink_info[, .(phenotype,  UniProt, HGNC.symbol, chr, gene_start, Panel, UKBPPP_ProteinID)],
       sep="\t", quote=FALSE, file="/hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/ml_inputs/v3/phenotype_info.txt")


olink_pheno = olink_pheno[phenotype %chin% unique(olink_SS$phenotype)]

#change again from : to _ in order to be easier to read in the next script 
olink_pheno$phenotype = gsub(":","_", olink_pheno$phenotype)
# write out 
fwrite(olink_pheno, sep="\t", quote=FALSE, file="/hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/ml_inputs/v3/phenotypes.txt")
