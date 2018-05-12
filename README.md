Rephine

RePhine is a method for identification of drug response related transcriptional regulators by modeling correlation patterns of targets. RePhine takes in pharmacogenomics data containing gene expression, copy number, mutation status, cancer type information and pharmacological profiles for estimation of associations between genes expression and drug response. It also takes in ChIP-seq data for target inference.

The RePhine method consists of the following steps:

Estimate the adjusted expression for each gene by removing the effects of copy number through linear regression.
Identify mutated genes related to drug response based on the CCLE mutation data and pharmacology profiles through an adaptive lasso model and likelihood ratio test.
Compute the partial correlation between a specific gene and a drug by controlling the confounding factors.
Repeat step 3 for next gene.
Down-weight non-cancer related genes.
Use likelihood ratio test to evaluate significance (uniP) of the TRs whose targets have concordant positive or negative partial correlations with drug response and apply the elastic net analysis combined with likelihood ratio test (multiP) to select the independent TRs.
Use a modified GSEA algorithm to further filter the top candidate lists derived from step 6. Regulatory potential scores across 1200 ChIP-seq data from Encode that are used in RePhine analysis are available on google drive.
Prerequisites

Python
Usage

First you should download dataset with:

$wget 
Results

Contact

liushiyi1991@hotmail.com
