import pandas as pd
import numpy as np
import joblib
import sys,os

def read_proteins_list(proteinomics_list_file):
    df = pd.read_csv(proteinomics_list_file,sep='\t')
    return list(df['phenotype'])


def read_protein_phenos(proteinomics_phenos_file, protein_name, sample_ids):
    df_pheno = pd.read_csv(proteinomics_phenos_file, delimiter='\t')
    df_pheno = df_pheno.loc[df_pheno['phenotype'] == protein_name]
    df_pheno = df_pheno.set_index('IID')
    return np.array(df_pheno.loc[sample_ids, 'value'])


def read_protein_genotypes(geno_file):
    df = pd.read_csv(geno_file, compression='gzip', sep='\t')
    sample_ids = list(df['varID'])
    df = df.set_index('varID')
    var_ids = list(df.columns)
    return sample_ids, var_ids


def run_experiments_5_folders_one_protein(proteinomics_list_file,proteinomics_genotype_path,results_path,models_path,protein_index_start,protein_index_end):

    #read the full list of protein unique ids & read the current protien id
    proteins_list = read_proteins_list(proteinomics_list_file)

    for protein_index in range(protein_index_start,protein_index_end):

        protein_name = proteins_list[protein_index]
        print("Start processing {}-{}".format(protein_index,protein_name))

        #read genotype matrix X and all sample ids and variants ids for the use of creating the model file (getting the corresponding varinat ID of extracted effect sizes)
        geno_file = proteinomics_genotype_path + protein_name + "_dosages.txt.gz"
        sample_ids,var_ids = read_protein_genotypes(geno_file)

        # Note folder 6 corresponding to the models trained with the whole PPP set
        for folder_count in range(1,7):

            model_file = models_path + protein_name + "_BR_model_" + str(folder_count) + ".pkl"

            br_model = joblib.load(model_file)
            effect_sizes = br_model.coef_
            data = {'rsid': var_ids, 'br_effect': effect_sizes}
            df = pd.DataFrame(data)

            effect_size_file = proteinomics_genotype_path + protein_name + "_variant_effects.txt"
            df1 = pd.read_csv(effect_size_file, sep='\t')

            df2 = pd.merge(df1, df, on='rsid')
            df3 = df2[["rsid", "chr", "pos", "effect_allele", "other_allele", "br_effect"]]

            results_file = results_path + protein_name + "_br_effects_" + str(folder_count) + ".txt"
            df3.to_csv(results_file, index=False, sep='\t')



if __name__ == "__main__":

        # The set of proteins to be processed in the protein list file in one job, i.e. the start and end index
        protein_index_start = int(sys.argv[1])
        protein_index_end = int(sys.argv[2])

        #The file lists all the protein IDs
        proteinomics_list_file = "/hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/ml_inputs/v3/phenotype_info_final.txt"

        #Folder of the genotype data (.txt file for genotype data matrix) as well as the GWAS based model files (_variant_effects.txt) generated previously
        proteinomics_genotype_path = "/hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/ml_inputs/v3/"

        #folder saves the .pkl models
        models_path = "/hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/ml_models/final/"

        # Folder saves the final model file (from .pkl filrs)
        results_path = "/hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/ml_models/final/O3K_models/"

        if os.path.isdir(results_path) == False:
            os.mkdir(results_path)

        run_experiments_5_folders_one_protein(proteinomics_list_file,proteinomics_genotype_path,results_path,models_path,protein_index_start,protein_index_end)

