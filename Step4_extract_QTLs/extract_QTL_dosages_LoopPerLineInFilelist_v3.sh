#!/bin/bash
#BSUB -q medium
#BSUB -n 5
#BSUB -J -extract_dosages_all_v3-
#BSUB -R "span[ptile=4]"
#BSUB -M 50000
#BSUB -o /hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/logs/extract_QTLs/extract_dosages_all_v3.j.o
#BSUB -e /hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/logs/extract_QTLs/extract_dosages_all_v3.j.e

out_dir=/hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/ml_inputs/v3
ls $out_dir/*_variant_effects.txt | sed "s#$out_dir/##" | sed 's/_variant_effects.txt//' > $out_dir/filelist.txt
input=$out_dir/filelist.txt
while IFS= read -r line;
        do
echo "$line"
phen="$line"
  # Determine the samples we need to keep for this phenotype
  grep "${phen}" $out_dir/phenotypes.txt | cut -f 1 > $out_dir/${phen}_IID.txt
  echo "#FID"$'\t'"IID" > $out_dir/${phen}.samples
  paste $out_dir/${phen}_IID.txt $out_dir/${phen}_IID.txt >> $out_dir/${phen}.samples
  rm $out_dir/${phen}_IID.txt

    # Get the list of pQTLs to extract
cut -f 1 $out_dir/${phen}_variant_effects.txt > $out_dir/${phen}_variant_ids.txt
    # Get their effect alleles
    grep -f $out_dir/${phen}_variant_ids.txt -w $out_dir/${phen}_variant_effects.txt | cut -f 1,4 > $out_dir/${phen}_variant_effect_alleles.txt 
    # Extract the dosages of the effect alleles
   /hpc/grid/wip_drm_targetsciences/users/dimitriou/Protein-PRS/test/plink2 \
--pfile /hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/UKB_500k_imputed_Broad_EUR_x_Olink_dedup_unambig_SNPs_maf0.005_ldthin0.8 \
           --out $out_dir/${phen}_dosages \
           --keep $out_dir/${phen}.samples \
           --extract $out_dir/${phen}_variant_ids.txt \
           --export A --export-allele $out_dir/${phen}_variant_effect_alleles.txt \
           --memory 50000 \
           --threads 5 \
           --silent

# remove temporary files
rm $out_dir/${phen}_variant_ids.txt 
rm $out_dir/${phen}_variant_effect_alleles.txt
rm $out_dir/${phen}_dosages.log

# change dosages.raw to dosages.txt
paste $out_dir/${phen}_dosages.raw > $out_dir/${phen}_dosages.txt 
rm $out_dir/${phen}_dosages.raw

# get the varID column and remove extra sample identifier information from each chromosome
module load ib
module load R/4.1.2
Rscript /hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/test_reformat_dosages_v3.R $out_dir $phen
rm $out_dir/${phen}_dosages.txt

# sample file no longer needed
rm $out_dir/${phen}.samples

done < "$input"

