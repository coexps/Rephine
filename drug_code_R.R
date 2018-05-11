read.table("CCLE_ACT.txt",sep="\t",header=T)->CCLE_drug #load CCLE drug activity from supp files
(CCLE_drug)[,c(1,3,13)]->CCLE_ACT #pick useful columns

####format transformation from cell line + drug +activity to matrix format
matrix(0,nrow=length(unique(CCLE_ACT[,1])),ncol=length(unique(CCLE_ACT[,2])))->drug_ACT
rownames(drug_ACT)<-unique(CCLE_ACT[,1])
colnames(drug_ACT)<-unique(CCLE_ACT[,2])
apply(as.matrix(CCLE_ACT),1,function(x) drug_ACT[x[1],x[2]]=x[3])
for(i in 1:length(CCLE_ACT[,1]))
{
    drug_ACT[as.vector(CCLE_ACT[i,1]),as.vector(CCLE_ACT[i,2])]<-as.numeric(as.vector(CCLE_ACT[i,3]))
}

####load CNV  and expression data of CCLE
read.table("ccle_cnv.txt",header=T,row.names=1)->CNVV
read.table("ccle_exp.txt",header=T,row.names=1)->EXPP

####pick the overlap samples
colnames(CNVV)->CNVV.colname
colnames(EXPP)->EXPP.colname
intersect(CNVV.colname, EXPP.colname)->inter.colname
rownames(CNVV)->CNVV.rowname
rownames(EXPP)->EXPP.rowname
intersect(CNVV.rowname, EXPP.rowname)->inter.rowname
rownames(drug_ACT)->cell_name
intersect(cell_name,inter.colname)-> inter.colname2
EXPP[inter.rowname, inter.colname2]-> EXPP2
CNVV[inter.rowname, inter.colname2]-> CNVV2

#####load mutation samples
read.delim("CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf",sep="\t")->maf
maf[grep("p.",maf[,"Protein_Change"]),]->maf2 #pick protein changing mutations
maf2[,c("Hugo_Symbol","Tumor_Sample_Barcode")]->maf3 # pick informative columns
table(maf3)->mut_pre ###transformate gene + sample into matrix
mut_pre[mut_pre>1]=1 # if one genes have two mutations in one sample, record only once.
mut_pre->mut
intersect(inter.colname2,colnames(mut))->inter.colname_mut
mut[,inter.colname_mut]->mut_mut
CNVV2[,inter.colname_mut]->CNVV_mut
drug_ACT[inter.colname_mut,]-> drug_ACT_mut
mut_mut[rowMeans(mut_mut)>0.05,]->mut_filter ###gene with too few mutation counts were removed

####adaptive lasso to select drug response related mutations that are independent of TR regulation
library(glmnet)
require(MASS)
#j="Erlotinib"

coef_all=c()
for(j in colnames(drug_ACT_mut))
{
    # format transformation
    rbind(mut_filter,unlist(drug_ACT_mut[,j]))->new_mut_mat
    t(new_mut_mat)->tnew_mut_mat
    na.omit(tnew_mut_mat)->tnew_mut_mat2
    
    X_var=tnew_mut_mat2[,-dim(tnew_mut_mat2)[2]] ##recurrent mutations
    y_var=tnew_mut_mat2[,dim(tnew_mut_mat2)[2]] ##response
    
    # calculate the weight of the adaptive lasso
    betals = ginv(t(X_var)%*% X_var)%*%t(X_var)%*%y_var
    w_alasso <- 1/abs(betals)
    
    # adaptive lasso to select mutations and use cross validation to determine the lambda
    cv.glmnet(X_var, y_var,penalty.factor = w_alasso,nfolds=dim(tnew_mut_mat2)[1],grouped=F)->res2
    
    #res2$lambda.min->bst_lam
    #res2$lambda.1se->bst_lam
    
    as.matrix(coef(res2, s = "lambda.min"))->coef1
    
    print(length(coef1))
    coef_all=cbind(coef_all,as.matrix(coef(res2, s = "lambda.min")))
    
    #plot(res2)
    #glmnet(X_var, y_var,family = "gaussian", alpha = 1, lambda= bst_lam,penalty.factor = w_alasso)->res_new
}

colnames(coef_all)<-colnames(drug_ACT_mut)

### use likelihood ratio test to evaluate the significance of univariate and multivariate association of selected mutation

feature_out=c()
for(jj in colnames(drug_ACT_mut))
{
    # format transformation
    rbind(mut_filter,unlist(drug_ACT_mut[,jj]))->new_mut_mat_sf
    t(new_mut_mat_sf)->tnew_mut_mat_sf
    na.omit(tnew_mut_mat_sf)->tnew_mut_mat_sf2
    
    lrt_res=c()
    setdiff(rownames(coef_all)[which(coef_all[,jj]!=0)],"(Intercept)")->pick_predictor #remove the Intercept columns
    # likelihood ratio test univariate and leave-one-out multivairate effect.
    if (length(pick_predictor)==0) {lrt_res=c()} else
    {
        for (tt in 1:length(pick_predictor)) #select one mutant in turns
        {
            fm1<-glm(tnew_mut_mat_sf2[,dim(tnew_mut_mat_sf2)[2]]~tnew_mut_mat_sf2[, pick_predictor]) #model of all selected mutants
            if(length(pick_predictor)==1) {fm2=glm(tnew_mut_mat_sf2[,dim(tnew_mut_mat_sf2)[2]]~NULL)} else
            {fm2<-glm(tnew_mut_mat_sf2[,dim(tnew_mut_mat_sf2)[2]]~tnew_mut_mat_sf2[, pick_predictor[-tt]])} #model of leave one mutant out
            fm3<-glm(tnew_mut_mat_sf2[,dim(tnew_mut_mat_sf2)[2]]~tnew_mut_mat_sf2[, pick_predictor[tt]]) # model of univariate
            lrtest(fm2,fm1)->loo  #Lrt of model of all selected mutants vs model of leave one mutant out
            lrtest(fm3)->single # lrt model of univariate
            c(loo[1,2]-loo[2,2],loo[2,5],single[1,2]-single[2,2],single[2,5],single[2,3])->lrt_1
            rbind(lrt_res, lrt_1)-> lrt_res
        }
        rownames(lrt_res)<-pick_predictor
    }
    c(feature_out, list(lrt_res))-> feature_out
}

# filter drug response related mutations based on the P-value of the likelihood ratio test
feature_filter_out=list()
for (i in 1:length(colnames(drug_ACT_mut)))

{
    # set p=0.01 as the cutoff
    (rownames(as.data.frame(feature_out[[i]])[which(feature_out[[i]][,1]<0&feature_out[[i]][,2]<0.01&feature_out[[i]][,3]>0&feature_out[[i]][,4]<0.01),]))->out1
    c(feature_filter_out, list(i=out1))-> feature_filter_out
    
}

# data transformation to match sample and drug with the results
feature_filter_out=list()
for (i in 1:length(colnames(drug_ACT_mut)))
{
    (rownames(as.data.frame(feature_out[[i]])[which(feature_out[[i]][,1]<0&feature_out[[i]][,2]<0.001&feature_out[[i]][,3]>0&feature_out[[i]][,4]<0.001),]))->out1
    c(feature_filter_out, list(i=out1))-> feature_filter_out
}
as.vector(unlist(feature_filter_out))-> feature_filter_out_all
matrix(0,ncol=length(colnames(drug_ACT_mut)),nrow=length(feature_filter_out_all))->feature_matrix
rownames(feature_matrix)<-feature_filter_out_all
colnames(feature_matrix)<-colnames(drug_ACT_mut)
for (i in 1:length(colnames(drug_ACT_mut)))
{
    feature_matrix [feature_filter_out[[i]],i]=1
}

# calculate the adjust expression by removing the confounding effect of CNVs through linear regression model
fun<-function(x)
{
    lm(unlist(EXPP2[x, inter.colname2])~unlist(CNVV2[x, inter.colname2]))[[1]][2]
}
apply(as.matrix((inter.rowname)),1,  fun)->CNV_lm

names(CNV_lm)<-inter.rowname
CNVV2* CNV_lm->CNV_factor
EXPP2-CNV_factor-> reg_factor


# calculate the partial correlation between adjusted expression and drug response by accounting for the confoundings.
sapply(inter.colname_mut,function(x) paste(unlist(strsplit(x,"_"))[-1],collapse="_"))->tissue #pick the tissue information of each cell line

tissue_factor=rep(1,length(tissue))
tissue_factor[which(tissue=="HAEMATOPOIETIC_AND_LYMPHOID_TISSUE")]=0 # diff blood cancer and solid tumors
names(tissue_factor)<-inter.colname_mut


fun2<-function(x) summary(glm(unlist(drug_ACT_mut[,j])~t(rbind(reg_factor[x,(inter.colname_mut)],mut_mut[rownames(feature_matrix)[feature_matrix[,j]!=0], inter.colname_mut], tissue_factor[inter.colname_mut]))
))$coefficients[2,3]
fun22<-function(x) {
    cbind(unlist(drug_ACT_mut[,j]),t(rbind(reg_factor[x,(inter.colname_mut)],mut_mut[rownames(feature_matrix)[feature_matrix[,j]!=0], inter.colname_mut], tissue_factor[inter.colname_mut])))->matrix_partial_cor
    na.omit(matrix_partial_cor)-> matrix_partial_cor
    pcor.test(matrix_partial_cor[,1],matrix_partial_cor[,2],matrix_partial_cor[,3:dim(matrix_partial_cor)[2]])->partial_cor
    (partial_cor[1:2])
} # define function for partial correlation coefficient calculation by pcor.test, var1 is the response, var2 is the adjusted expression, var3 is the tumor type, var4 is the selected mutation information

# data transformation from list to matrix.
stage_all_partial=c()
partial_P=c()
for (j in colnames(drug_ACT_mut))
{
    print(j)
    sapply(rownames(reg_factor), fun22)->partial_output
    
    cbind(stage_all_partial, partial_output[1,])-> stage_all_partial
    cbind(partial_P, partial_output[2,])-> partial_P
    
}
stage_all_partial2<-matrix(0,nrow=dim(stage_all_partial)[1],ncol=dim(stage_all_partial)[2])
for(i in 1:dim(stage_all_partial)[2])
{
    unlist(stage_all_partial[,i])->stage_all_partial2[,i]
}
colnames(stage_all_partial2)<-colnames(drug_ACT_mut)
rownames(stage_all_partial2)<-rownames(stage_all_partial)



##### integre TCGA data by removing non-differential genes between cancer and tumor across TCGA cancer types
read.table("merge",sep="\t",header=T,row.names=1)->diff
gsub(" ","",rownames(diff))->diff2_name
rownames(diff)<-diff2_name
apply(diff,1,function(x) sort(abs(x[-15]),)[5]<=5e-3)->tttt # if one gene is differential less than 1/3 cancer types, set weight to 0 only a few genes were filtered in this step but it is quite useful to remove the random bias
as.numeric(tttt)->tttt_num
names(tttt_num)<-names(tttt)

apply((stage_all_partial2)[inter.name_chip,],2,function(x) x* tttt_num[inter.name_chip])->stage_all_tmp
stage_all_tmp-> stage_all_partial2


# ChIP-seq data selections from the replicates

read.table("CHIP_SEQ/merge.txt",header=T)->merge_score #load RP scores across different ChIP-seq data
intersect(unlist(merge_score[,2]),rownames(reg_factor))->inter.name_chip
sort(inter.name_chip)->sort_inter.name_chip
single_TF_filter=c()
sort(inter.name_chip)->sort_inter.name_chip
# format transformation: 1)remove the duplicate transcripts, 2) match genes in expression data to those in ChIP-seq data
for (n in out_names)
{
    
    order(merge_score[,n],decreasing=T)->order_merge
    merge_score[order_merge,2]->order_genes
    as.vector(order_genes)-> order_genes
    which(!duplicated(order_genes))->uniq_gene_order
    merge_score[order_merge [uniq_gene_order],n]->uniq_value
    names(uniq_value)<-order_genes[uniq_gene_order]
    uniq_value[sort_inter.name_chip]-> uniq_value2
    cbind(single_TF_filter, uniq_value2)-> single_TF_filter
}

colnames(single_TF_filter)<-out_names
stage_all_partial2[sort_inter.name_chip,drug_name]->yy
#sapply(out_names,function(x) lrtest(glm(yy~ single_TF_filter[,x]))[2,5])->TF_filter_P

sapply(out_names,function(x) {glm(yy~ single_TF_filter[,x])->func_;c(func_$coefficients[2],lrtest(func_)[2,5])})->TF_coef_out
t(TF_coef_out)-> TF_coef_out2 #calculate the significance of each replicate and also the univariate significance.
TF_coef_out2[,2]-> TF_filter_P
names(TF_filter_P)<- out_names
sort(TF_filter_P)->sort_TF_filter_P




split_TF<-function(TF_names)
{
    strsplit(TF_names[1],"\\.")  ->split_1
    split_1[[1]][length(unlist(split_1))-2]-> split_2
    strsplit(split_2,"_")->split_2_2
    split_2_2[[1]][length(unlist(split_2_2))]  ->split_3
    
} #define function to pick the TR names from the files
sapply(names(sort_TF_filter_P), split_TF)-> sort_TF_filter_P_pick
cbind(sort_TF_filter_P, sort_TF_filter_P_pick)->merge_sort_TF_filter_P
merge_sort_TF_filter_P [which(!duplicated(merge_sort_TF_filter_P[,2])),]-> merge_sort_TF_filter_P_filter #select the most significant replicate for the following analysis
merge_sort_TF_filter_P_filter[!merge_sort_TF_filter_P_filter[,2]%in%c("CTCF","POLR2A","POLR2AphosphoS5","POLR2AphosphoS2"),]-> merge_sort_TF_filter_P_filter2 #remove CTCF POLR2 factors into consideration because they are not the traditional TRs

rownames(merge_sort_TF_filter_P_filter2)-> TF_chipseq


single_TF_filter[,TF_chipseq]-> merge_TF_pick
colnames(merge_TF_pick)->filter_feature_TFs
as.vector(filter_feature_TFs)-> filter_feature_TFs
merge_TF_pick->X_var_TF
stage_all_partial2[sort_inter.name_chip,drug_name]->Y_var_TF


##### permutation function
# define the function to calculate the enrichment scores of the modified GSEA
fun3<-function(mm, Par_Score, RP_score)
{
    RP_score_adj<-RP_score/max(RP_score)
    sum(abs(Par_Score[inter.name_chip]^mm)* RP_score_adj[inter.name_chip])->NR_sum
    sort(Par_Score,decreasing=T)->sort_Par_Score
    names(sort_Par_Score)->sort_name
    RP_score_adj[sort_name]-> sort_RP_score_adj
    P_hit<-sort_RP_score_adj*(abs(sort_Par_Score)^mm)/NR_sum
    Non_Sum<-length(RP_score_adj)-sum(sort_RP_score_adj)
    #P_miss<-(1-P_gene)/(1-mean(P_Chip))/length(P_Chip)
    P_miss<-(1-sort_RP_score_adj)/Non_Sum
    #P_miss
    #P_hit
    P_hit-P_miss
}
# define the permutation function to calculate the enrichment scores of the modified GSEA in random cases
fun_random<-function(mm, Par_Score, RP_score)
{
    RP_score_adj<-RP_score/max(RP_score)
    sort(Par_Score,decreasing=T)->sort_Par_Score
    names(sort_Par_Score)->sort_name
    #RP_score_adj[sort_name]-> sort_RP_score_adj
    sample(RP_score_adj)-> sort_RP_score_adj
    sum(abs(sort_Par_Score^mm)* sort_RP_score_adj)->NR_sum
    P_hit<-sort_RP_score_adj*(abs(sort_Par_Score)^mm)/NR_sum
    Non_Sum<-length(RP_score_adj)-sum(sort_RP_score_adj)
    #P_miss<-(1-P_gene)/(1-mean(P_Chip))/length(P_Chip)
    P_miss<-(1-sort_RP_score_adj)/Non_Sum
    #P_miss
    #P_hit
    P_hit-P_miss
}


## elastic net for TR selections

cv.glmnet(X_var_TF, Y_var_TF, alpha=0.8, nfolds=160 )-> elasticnet #we set alpha fixed to 0.8 to save the calculation time
as.matrix(coef(elasticnet,s = "lambda.1se"))->elastic_coef # pick 1se lambda to avoid over-fit
rownames(elastic_coef)[which(elastic_coef!=0)]-> elastic_coef_names
times=1000 # set permutation times
fmTF1<-glm(Y_var_TF ~ X_var_TF[, elastic_coef_names[-1]]) # model of all selected TRs
lrt_TFres=c()
# likelihood ratio test & permutation test
for (tt in elastic_coef_names[-1])
    {
        # likelihood ratio test
        setdiff(elastic_coef_names[-1],tt)-> elastic_coef_names_out #leave one TR out in turn
        print(tt)
        #fmTF1<-glm(Y_var_TF ~ X_var_TF[, elastic_coef_names[-1]])
        fmTF2<-glm(Y_var_TF ~ X_var_TF[, elastic_coef_names_out ]) # model of TRs with leaving one out
        #fmTF3<-glm(Y_var_TF ~ X_var_TF[, tt ])
        #lrtest(fmTF3)->single
        lrtest(fmTF1, fmTF2)->loo #likelihood ratie test of model of all selected TRs &  TRs with leaving one out to detect the multivariate effect.
        #c(fmTF3$coefficients[2],loo[1,2]-loo[2,2],loo[2,5],single[2,4],single[2,5])->lrt_TF1
        rand_max=c()
        rand_min=c()
        # permutation test
        for(i in 1: times)
        {
        cumsum(fun_random(1,Y_var_TF,X_var_TF[, tt]))->cum_rand
        rand_max=c(rand_max,max(cum_rand))
        rand_min =c(rand_min,min(cum_rand))
        }
        cumsum(fun3(1,Y_var_TF,X_var_TF[, tt]))->cum_value
    
        value_max=max(cum_value)
        value_min = min(cum_value)
        Perm_P=min(length(which(value_max< rand_max))/times,length(which(value_min >  rand_min))/times)
        c(loo[2,4],loo[2,5], Perm_P)-> lrt_TF1
        rbind(lrt_TFres, lrt_TF1)-> lrt_TFres
    }
rownames(lrt_TFres)<-elastic_coef_names[-1]
#lrt_TFres[which(lrt_TFres[,3]<0.001& lrt_TFres[,5]<0.00001),]-> merge_TFs_filter
TF_coef_out2 [TF_chipseq,]->output_1
cbind(output_1[TF_chipseq,], elastic_coef[TF_chipseq,])->output_2 #combine loo significance and univariate signaficance
merge(output_2, lrt_TFres,by.x=0,by.y=0,all.x=T)->out_last
#out_last[which(out_last[,6]<0.01& out_last[,3]<1e-3),]-> filter_out_last
#out_last[which(out_last[,6]<1& out_last[,3]<1),]-> filter_out_last2
write.table(out_last,paste("CCLE_drug/",drug_name,".txt",sep=""),quote=F,sep="\t") # output the results of each drugs

}




















