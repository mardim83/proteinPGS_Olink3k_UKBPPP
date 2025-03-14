# proteinPGS_Olink3k_UKBPPP

This repository deposits codes used to train genetic scores of proteomics Olink3k using UK-Biobank Pharma Proteomics Project (UKB-PPP) that utilized in the study: Dimitriou et al. Protein-Polygenic scores to identify etiological risk factors in cancer. Manuspript in preparation (2025).
For training of the protein-PGS models we utilized published proteomics data for 2940 proteins as well as GWAS results for the same proteins from the UKB-PPP discovery study of 34,557 participants of European ancestry as previously described (Sun et al., Nature 2023). We selected the Baysian Ridge (BR) methodology with the previously optimized variant-selection scheme of P < 5×10−8 on genome-wide variants as previously described (Xu et al., Nature 2023).
Briefly, prior to the training we initially performed appropriate variant selection by utilizing the imputed and QCed genotypes. Dublicated variants removed, variants in multi-allelic sites removed, strand ambiguous alleles removed, variants LD thinned with MAF > 0.5% at R2=0.8 using the indep-pairwise method implemented in plink2 and finally kept only variants with pQTLs with P < 5×10−8 from their protein-GWAS result. By utilizing the summary statistics from GWAS for all the available Olink proteins for the UKB discovery study and using only variants with above-mentioned P value threshold we were able to train models for 2601 traits out of the 2940 traits (individual OlinkIDs) with available proteomic data. For the training of the models with the BR methodology the scikit-learn package was implemented and scipy package for validation analysis to derive the explained variance (R2) and Spearman’s rank correlation coefficient (Rho) with Python (v3.7.12). As internal validation we used a 5-fold cross validation setting to equally partition the UKB discovery study in 5 subsets (80% samples for training – 20% samples for validation) which resulted in five different genetic-score models per trait. To confirm superiority of BR over the P+T in our UKB training dataset we use their R2 as performance metrics from all the possible genetic score from this partition. Then the mean of R2 in the 5 testing samples in UKB was utilized for subsequent performance comparisons with the external validation study in Finngen. Given the overall consisted performance between the 5 cross-validated models per trait we finally produced one model per trait (2601 models) using the full UKB discovery study of the 34,557 participants of European ancestry. To validate these models externally we applied them in 1962 individuals from Finngen with available proteomics data from the Olink3k explore panel.

# Folder contents description
*	Genetic score development for Olink3k proteomics traits:
    * Step1_dedup: Annotate variants with unique identifiers to flag and remove duplicates;
    * Step2_ldthin: remove multi-allelic, ambiguous (A/T, G/G) variants and variants with a MAF < 0.5%, and ld-thin variants with r2=0.8 (i.e. indep-pairwise 1000kb 0.8 in plink2);
    * Step3_collate_QTLs: collate the list of QTLs for variant selection from GWAS summary statistics;
    * Step4_extract_QTLs: select the list QTLs with p<5×10−8, and extract their dosages of the effect alleles as input data of Bayesian Ridge;
    * Step5_genetic_score_training: training genetic score models using Bayesian ridge.
