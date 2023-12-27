# MI-SRT-Correlation

In this study, we apply single imputation (SI) and multiple imputation (MI) on twelve pairs of scRNA-seq and SRT datasets. We use the data curated by Li et al. (2022) in the Google Drive repo: https://drive.google.com/drive/folders/1pHmE9cg_tMcouV1LFJFtbyBJNp7oQo9J?usp=sharing. In the Google Drive repo, the imaging-based SRT data can be found in DataUpload.zip after extraction. The twelve pairs of SRT and scRNA-seq datasets we use are in subfolders named as Dataset1, 2, 4, 5, 7, 10, 11, 12 14, 15, 16, and 17. 

In this Github repo, we present the scripts used for imputation experiments and performance comparison for the pair of SRT and scRNA-seq datasets in Dataset1, where the SRT dataset is named "Insitu_count.txt" and the scRNA-seq dataset is named "scRNA_count.txt". The SRT dataset is measured using the seqFISH technology while the scRNA-seq dataset is measured using the 10X Chromium technology. Therefore, we named this pair as "seqFISH_10X".

Predict_gene.txt: We randomly select 100 sets of 35 out of the 347 genes measured in the SRT dataset and store 100 lines of corresponding gene names in this txt file.

seqFISH_10X.py: This python script takes the Predict_gene.txt, a scRNA_count.txt and a Insitu_count.txt as input. For single imputation, the script generates "Fold_[k]_SI.txt" for each randomly selected set of 35 genes, where k=1,..,100. For multiple imputation, the script generates "Fold_[i]_MI_[j].txt" for each set of genes and 100 random imputations, where i=1,...,100 and j=1,...,100.

ComputeCorrelation.R: This R script generates correlation matrices for original data, single imputation data, and multiple imputation data in pearson, spearman, and kendall correlation. The files are named as "Fold_[k]_Original_[cor.method].txt", "Fold_[k]_SI_[cor.method].txt", and "Fold_[i]_MI_[j][cor.method].txt"
 