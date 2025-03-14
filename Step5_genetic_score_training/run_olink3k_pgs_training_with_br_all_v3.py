import pandas as pd
import numpy as np
import datetime
from sklearn.model_selection import KFold
from BayesianRidge import full_fit_BayesianRidge
from sklearn.linear_model import BayesianRidge
from Traditional_GRS import traditional_GRS_selected_vars
import joblib
import sys
import os
import time


def read_proteins_list(proteinomics_list_file):
    df = pd.read_csv(proteinomics_list_file,sep='\t')
    return list(df['phenotype'])


def read_protein_phenos(proteinomics_phenos_file, protein_name, sample_ids):
    df_pheno = pd.read_csv(proteinomics_phenos_file, delimiter='\t')
    df_pheno = df_pheno.loc[df_pheno['phenotype'] == protein_name]
    df_pheno = df_pheno.set_index('IID')
    return np.array(df_pheno.loc[sample_ids, 'value'])


def read_protein_genotypes(geno_file):
    df = pd.read_csv(geno_file,compression='gzip', sep='\t')
    sample_ids = list(df['varID'])
    df = df.set_index('varID')
    var_ids = list(df.columns)
    X = np.array(df.loc[:,:])
    return sample_ids,var_ids,X

def full_fit_BayesianRidge_no_test(X, y, alpha_1, alpha_2, lambda_1, lambda_2):
    model = BayesianRidge(alpha_1=alpha_1, alpha_2=alpha_2, lambda_1=lambda_1, lambda_2=lambda_2)
    model.fit(X, y)
    return model

def run_experiments_5_folders_one_protein(proteinomics_list_file,proteinomics_genotype_path, proteinomics_phenos_file,results_path,models_path,beta_path,protein_index, alpha_1, alpha_2, lambda_1, lambda_2):

    #read the full list of protein (or other type of omic trait) unique ids & read the current protein (or ther type of trait) id
    proteins_list = read_proteins_list(proteinomics_list_file)
    protein_name = proteins_list[protein_index]

    print("Start processing {}-{}".format(protein_index,protein_name))

    #read genotype matrix X and all sample ids and variants ids
    geno_file = proteinomics_genotype_path + protein_name + "_dosages.txt.gz"
    sample_ids,var_ids,X = read_protein_genotypes(geno_file)

    print("Number of Variants {}".format(len(var_ids)))
    print("Number of Samples {}".format(len(sample_ids)))

    #read protein levels of the given proteins and all samples
    y = read_protein_phenos(proteinomics_phenos_file, protein_name, sample_ids)

    results_file = results_path + protein_name + "_BR_UNI_prs.txt"

    f = open(results_file,'w')
    f.write('Time\tProtein\tN_Vars\tFolder\tBR_r2\tBR_sr\tUNI_r2\tUNI_sr\n')

    folder_count = 0
    kf = KFold(n_splits=5, shuffle=True, random_state=21)
    for train_index, test_index in kf.split(y):
        folder_count += 1

        print("folder-{}".format(folder_count))

        x_train, x_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]

        print("{}-folder-{} Running BayesianRidge...".format(datetime.datetime.now(), folder_count))
        BR_model,br_r,br_r2,br_envs,br_sr = full_fit_BayesianRidge(x_train, x_test, y_train, y_test, alpha_1, alpha_2, lambda_1, lambda_2)
        print("r: {}, r2: {}, env: {}, sr: {}".format(br_r, br_r2,br_envs,br_sr))
        model_file = models_path + protein_name + "_BR_model_" + str(folder_count) + ".pkl"
        joblib.dump(BR_model, model_file)

        print("{}-folder-{} Running Univariant method...".format(datetime.datetime.now(), folder_count))
        beta_file = beta_path + protein_name + "_variant_effects.txt"
        grs_r, grs_r2,grs_env,grs_sr = traditional_GRS_selected_vars(beta_file, x_test, y_test, var_ids)
        print("r: {}, r2: {}, env: {}, sr: {}".format(grs_r, grs_r2,grs_env, grs_sr))

        write_text = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(datetime.datetime.now(),protein_name,X.shape[1],folder_count,br_r**2,br_sr,grs_r**2, grs_sr)
        f.write(write_text)
        f.flush()
    
    folder_count=6
    print("{}-folder-{} Running BayesianRidge - full dataset...".format(datetime.datetime.now(), folder_count))
    BR_model = full_fit_BayesianRidge_no_test(X, y, alpha_1, alpha_2, lambda_1, lambda_2)
    model_file = models_path + protein_name + "_BR_model_" + str(folder_count) + ".pkl"
    joblib.dump(BR_model, model_file)
    write_text = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(datetime.datetime.now(),protein_name,X.shape[1],folder_count,float("nan"),float("nan"),float("nan"), float("nan"))
    f.write(write_text)
    f.flush()

    f.close()

if __name__ == "__main__":
    

        # the omic trait index in a platform
        protein_index = int(sys.argv[1])
        
        # BR priors for traits in the platform
        alpha_1 = float(sys.argv[2])
        alpha_2 = float(sys.argv[3])
        lambda_1= float(sys.argv[4])
        lambda_2 = float(sys.argv[5])
        
        
        # trait level file
        proteinomics_phenos_file = "/hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/ml_inputs/v3/phenotypes_final.txt"
        
        # trait list on the platform
        proteinomics_list_file = "/hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/ml_inputs/v3/phenotype_info_final.txt"
        
        #folder store all genotype files at each folder
        proteinomics_genotype_path = "/hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/ml_inputs/v3/"
        
        # results path
        results_path = "/hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/results/final/"

        if os.path.isdir(results_path) == False:
            os.mkdir(results_path)
        
        # model path
        models_path = "/hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/ml_models/final/"

        if os.path.isdir(models_path) == False:
            os.mkdir(models_path)
        
        # path stores the the betas of selected QLT variants from GWAS
        beta_path = "/hpc/grid/wip_drm_targetsciences/users/dimitriou/OMICSPRED/Olink_3k/ml_inputs/v3/"

        
        # model training for one omic trait
        run_experiments_5_folders_one_protein(proteinomics_list_file,proteinomics_genotype_path, proteinomics_phenos_file,results_path,models_path,beta_path,protein_index, alpha_1, alpha_2, lambda_1, lambda_2)

