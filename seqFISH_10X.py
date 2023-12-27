# Load packages
import scanpy as sc
import pandas as pd
import numpy as np
import copy
import pickle

import scvi
from scvi.model import GIMVI

# Prepare data
PATH = 'Your work directory'
RNA_path = PATH + '/Your scRNA-seq dataset.txt'
Spatial_path = PATH + '/Your SRT dataset.txt'
Predict_path = PATH + '/Predict_gene.txt'
Output_path = PATH

RNA_data_adata = sc.read(RNA_path, sep = '\t', first_column_names = True)
Spatial_data_adata = sc.read(Spatial_path, sep = '\t', first_column_names = False)

predict_gene_df = pd.read_table(Predict_path, header = None)
predict_gene = []
n_fold = 100
for i in range(n_fold):
    predict_list = list(predict_gene_df.iloc[i,:])
    predict_gene.append(predict_list)

# Define gimVI training function
def gimVI_Train(k, distributions, library):
    global RNA_data_adata, Spatial_data_adata, predict_gene
    
    test_list = predict_gene[k]
    Genes = list(Spatial_data_adata.var_names)
    test_gene_idx = [Genes.index(x) for x in test_list]
    test_genes = np.array(Genes)[test_gene_idx]
    n_genes = len(Genes)
    train_gene_idx = sorted(set(range(n_genes)) - set(test_gene_idx))
    train_genes = np.array(Genes)[train_gene_idx]
    
    seq_data = copy.deepcopy(RNA_data_adata)
    seq_data = seq_data[:, Genes].copy()
    sc.pp.filter_cells(seq_data, min_counts = 1)
    scvi.data.setup_anndata(seq_data)
    
    spatial_data = copy.deepcopy(Spatial_data_adata)
    spatial_data_partial = spatial_data[:, train_genes].copy()
    sc.pp.filter_cells(spatial_data_partial, min_counts = 1)    
    scvi.data.setup_anndata(spatial_data_partial)
    
    model = GIMVI(seq_data, spatial_data_partial, 
                  generative_distributions = distributions, 
                  model_library_size = library)
    
    model.train(n_epochs = 400)
    
    return model, test_gene_idx, spatial_data_partial.obs_names   

# Output single and multiple prediction
n_fold = 100
n_pred = 100

for k in range(n_pred):
    model, pred_gene_idx, cell_names = gimVI_Train(k, ['zinb', 'nb'], [True, False])
    _, imputation_si = model.get_imputed_values(deterministic = True, normalized = False)

    Genes = list(Spatial_data_adata.var_names)
    pred_genes = np.array(Genes)[pred_gene_idx]  
    si_result = pd.DataFrame(imputation_si[:, pred_gene_idx], columns = pred_genes)
    si_result.to_csv(Output_path + '/Fold_' + str(k+1) + '_SI.txt', header = True, index = False, sep = '\t')

    original_data = pd.DataFrame(Spatial_data_adata[cell_names].X.astype(int), columns = Genes)
    original_data.to_csv(Output_path + '/Fold_' + str(k+1) + '_Original.txt', header = True, index = False, sep = '\t')    

    for i in range(n_pred):
        _, imputation_mi = model.get_imputed_values(deterministic = False, normalized = False)
        mi_result = pd.DataFrame(imputation_mi[:, pred_gene_idx], columns = pred_genes)
        mi_result.to_csv(Output_path + '/Fold_' + str(k+1) + '_MI_' + str(i+1) + '.txt', header = True, index = False, sep = '\t')