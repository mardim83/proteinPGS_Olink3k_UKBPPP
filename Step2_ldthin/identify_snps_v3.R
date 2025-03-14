setwd("/hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k")
getwd()

library(data.table)

pvar = fread("/hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/UKB_500k_imputed_Broad_EUR_x_Olink_dedup.pvar")

#test what the sript does for the first 100000 snps (only in interactive sessions)
#pvar = head(pvar, n=100000L)

# Remove multi-allelic sites:
pvar$location = paste(pvar$`#CHROM`, pvar$POS, sep= (":") )
#sum(duplicated(pvar$location))

multi = pvar[grepl("*", location), .N, by=location][N > 1]
pvar = pvar[!multi, on = .(location)]


# Function for flipping the strand of an allele.
# Uses a series of gsub calls to replace A's with T's,
# G's with C's, and vice-versa. Also works for alleles
# with more than one nucleotide (e.g. indels).
flip_strand <- function(x) {
  # Swap each letter for a dummy, we need this intermediate
  # step so we can distinguish between alleles when swapping.
  # E.g if we did A -> T then T -> A we'd end up with all A's
  # and no T's. instead we do A -> V -> T and T -> X -> A.
  x <- gsub("A", "V", x)
  x <- gsub("T", "X", x)
  x <- gsub("C", "Y", x)
  x <- gsub("G", "Z", x)
  x <- gsub("V", "T", x)
  x <- gsub("X", "A", x)
  x <- gsub("Y", "G", x)
  x <- gsub("Z", "C", x)
  return(x)
}

# Remove strand ambiguous alleles:
#removed = pvar[REF == flip_strand(ALT)]
pvar = pvar[REF != flip_strand(ALT)]


# Filter to SNPs
pvar = pvar[nchar(REF) == 1 & nchar(ALT) == 1]

pvar$location = NULL

# Write out list of variants to extract prior to LD-thinning
fwrite(pvar[,.(ID)], col.names=FALSE, quote=FALSE, file="/hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/ldthinned/keep_v3.txt")
