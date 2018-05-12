# Rephine

RePhine is a method for identification of drug response related transcriptional regulators by modeling correlation patterns of targets. RePhine takes in pharmacogenomics data containing gene expression, copy number, mutation status, cancer type information and pharmacological profiles for estimation of associations between genes expression and drug response. It also takes in ChIP-seq data for target inference.

The RePhine method consists of the following steps:
1)	Estimate the adjusted expression for each gene by removing the effects of copy number through linear regression.
2)	Identify mutated genes related to drug response based on the CCLE mutation data and pharmacology profiles through an adaptive lasso model and likelihood ratio test.
3)	Compute the partial correlation between a specific gene and a drug by controlling the confounding factors.
4)	Repeat step 3 for next gene.
5)	Down-weight non-cancer related genes.
6)	Use likelihood ratio test to evaluate significance (uniP) of the TRs whose targets have concordant positive or negative partial correlations with drug response and apply the elastic net analysis combined with likelihood ratio test (multiP) to select the independent TRs.
7)	Use a modified GSEA algorithm to further filter the top candidate lists derived from step 6.


## Prerequisites
* Python
* R and packages "glmnet", "MASS", "lmtest", "ppcor", "psych", "ROCR"

## Documents
RP_calculation_modified.py is the modified python script which calculates the RP scores by consider both distance and peak strength. It takes in peak bed files of ChIP-seq data.

permutation_code.R is the simulation R code. it simulate TR and targets expression data, noise data and RP scores.

drug_code_R.R is the RePhine R code. 

## Resources
Regulatory potential scores across 1200 ChIP-seq data from Encode that are used in RePhine analysis are available on [google drive](https://drive.google.com/open?id=1wWCKDMSbTBpfuw_s2pYf4RPw7QqjDDIr).

## Contact
_wangxujun87@sjtu.edu.cn_

