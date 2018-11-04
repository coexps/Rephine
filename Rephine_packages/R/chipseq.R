chipseq <- function(chip=matrix(0,0,0),stage_all_partial=matrix(0,0,0), genes=NULL, samples=NULL, ...)
{
  intersect(unlist(chip[,2]),rownames(stage_all_partial))->inter.name_chip
  sort(inter.name_chip)->sort_inter.name_chip
  # print(sort_inter.name_chip)
  apply(chip[,3:(dim(chip)[2])],2,function(x) {sum(x==0)/length(x)})->qc_ratio #filter too many peaks with noise
  # print (qc_ratio)
  names(which(qc_ratio>0.1))->out_names
  #print (out_names)
# Check reg_factor and input matrix
#single_TF_filter<-apply(x[,out_names],2, function(x) .uniq_v(x,genes))
    single_TF_filter=c()
    for (n in out_names)
    {
        #print (n)
        order(chip[,n],decreasing=T)->order_merge
        chip[order_merge,2]->order_genes
        as.vector(order_genes)-> order_genes
        which(!duplicated(order_genes))->uniq_gene_order
        chip[order_merge [uniq_gene_order],n]->uniq_value
        names(uniq_value)<-order_genes[uniq_gene_order]
        uniq_value[sort_inter.name_chip]-> uniq_value2
        #print(uniq_value2)
        cbind(single_TF_filter, uniq_value2)-> single_TF_filter
    }


    colnames(single_TF_filter)<-out_names
    out_last2=c()
    #print (colnames(stage_all_partial))
    for (drug_names in colnames(stage_all_partial))

    {
        print (drug_names)
          stage_all_partial[sort_inter.name_chip,drug_names]->yy
          sapply(out_names,function(x) {glm(yy~ single_TF_filter[,x])->func_;c(func_$coefficients[2],lrtest(func_)[2,5])})->TF_coef_out
          t(TF_coef_out)-> TF_coef_out2 #calculate the significance of each replicate and also the univariate significance.
          TF_coef_out2[,2]-> TF_filter_P
          names(TF_filter_P)<- out_names
          sort(TF_filter_P)->sort_TF_filter_P

#   print (TF_coef_out)
          T_filter_P<-t(TF_coef_out)[,2]
          #print(names(sort_TF_filter_P))
          sort_TF_filter_P_pick<-sapply(names(sort_TF_filter_P), .splite_TF)
          cbind(sort_TF_filter_P, sort_TF_filter_P_pick)->merge_sort_TF_filter_P
          #print (merge_sort_TF_filter_P)
          merge_sort_TF_filter_P [which(!duplicated(merge_sort_TF_filter_P[,2])),]-> merge_sort_TF_filter_P_filter #select the most significant replicate for the following analysis
          # print (merge_sort_TF_filter_P_filter)
          merge_sort_TF_filter_P_filter[!merge_sort_TF_filter_P_filter[,2]%in%c("CTCF","POLR2A","POLR2AphosphoS5","POLR2AphosphoS2"),]-> merge_sort_TF_filter_P_filter2 #remove CTCF POLR2 factors into consideration because they are not the traditional TRs
            rownames(merge_sort_TF_filter_P_filter2)-> TF_chipseq
            single_TF_filter[,TF_chipseq]-> merge_TF_pick
            colnames(merge_TF_pick)->filter_feature_TFs
            as.vector(filter_feature_TFs)-> filter_feature_TFs
            merge_TF_pick->X_var_TF
            stage_all_partial[sort_inter.name_chip,drug_names]->Y_var_TF
            #print (Y_var_TF)
            #print (X_var_TF)

       cbind(Y_var_TF,X_var_TF)->merge_Var
              as.data.frame(merge_Var)->merge_Var2
       cv.glmnet(X_var_TF, Y_var_TF, alpha=0.8, nfolds=dim(X_var_TF)[2] )-> elasticnet
       as.matrix(coef(elasticnet,s = "lambda.1se"))->elastic_coef
       #print (elastic_coef)
       #print (elastic_coef)
       if (length(which(elastic_coef[-1]!=0)>0))
       {
           rownames(elastic_coef)[which(elastic_coef!=0)]-> elastic_coef_names
           times=1000
           fmTF1<-glm(Y_var_TF ~ X_var_TF[, elastic_coef_names[-1]])
           #print (fmTF1)
           lrt_TFres=c()
           for (tt in elastic_coef_names[-1])
           {
               setdiff(elastic_coef_names[-1],tt)-> elastic_coef_names_out
               # print (elastic_coef_names_out)
               # print(tt)
               #fmTF1<-glm(Y_var_TF ~ X_var_TF[, elastic_coef_names[-1]])
               if (length(elastic_coef_names_out)>0){fmTF2<-glm(Y_var_TF ~ X_var_TF[, elastic_coef_names_out ])}
               else {fmTF2<-glm(Y_var_TF ~NULL)}
               #fmTF3<-glm(Y_var_TF ~ X_var_TF[, tt ])
               #lrtest(fmTF3)->single
               lrtest(fmTF1, fmTF2)->loo
               #c(fmTF3$coefficients[2],loo[1,2]-loo[2,2],loo[2,5],single[2,4],single[2,5])->lrt_TF1
               rand_max=c()
               rand_min=c()
               for(i in 1: times)
               {
                   cumsum(.permutation_random(1,Y_var_TF,X_var_TF[, tt]))->cum_rand
                   rand_max=c(rand_max,max(cum_rand))
                   rand_min =c(rand_min,min(cum_rand))

               }
               cumsum(.permutation_modified(1,Y_var_TF,X_var_TF[, tt]))->cum_value

               value_max=max(cum_value)
               value_min = min(cum_value)
               Perm_P=min(length(which(value_max< rand_max))/times,length(which(value_min >  rand_min))/times)

               c(loo[2,4],loo[2,5], Perm_P)-> lrt_TF1
               rbind(lrt_TFres, lrt_TF1)-> lrt_TFres


           }
       rownames(lrt_TFres)<-elastic_coef_names[-1]
       #lrt_TFres[which(lrt_TFres[,3]<0.001& lrt_TFres[,5]<0.00001),]-> merge_TFs_filter


       TF_coef_out2 [TF_chipseq,]->output_1
       cbind(output_1[TF_chipseq,], elastic_coef[TF_chipseq,])->output_2

       merge(output_2, lrt_TFres,by.x=0,by.y=0,all.x=T)->out_last
       rownames(out_last)<-out_last[,1]
       out_last[,-1]->out_last
       #print (cbind(out_last,drug_names))
       #print (head(cbind(out_last,drug_names)))
       rbind(out_last2,as.matrix(cbind(out_last,drug_names)))->out_last2
       #print (out_last2)
       #print (head(out_last2))
       }
       else
       {

           TF_coef_out2 [TF_chipseq,]->output_1
           cbind(output_1[TF_chipseq,], elastic_coef[TF_chipseq,])->output_2

           cbind(output_2, NA,NA,NA)->out_last
           #print (out_last)
           #print (cbind(out_last,drug_names))
           #print (dim(cbind(out_last,drug_names)))
           #print(cbind(out_last,drug_names))
           #print (out_last2)
           rbind(out_last2,as.matrix(cbind(out_last,drug_names)))->out_last2
           #print (head(out_last2))

#print (out_last)
       }
    }
    out_last2
}
#.uniq_v<-function(x,y){
# genes<-y
# x<-x
# x_ord<-order(x,decreasing = TRUE)
# gene_ord<-as.vector(genes[x_ord])
# gene_uniq_ord<-which(!duplicated(gene_ord))
# value_uniq<-x[x_ord[gene_uniq_ord]]
# names(value_uniq)<-gene_ord[gene_uniq_ord]
# value_uniq
#}
.splite_TF<-function(x)
{
    split1<-strsplit(x[1],"\\.")
    split2<-split1[[1]][length(unlist(split1))-2]
    split2_2<-strsplit(split2,"_")
    split3<-split2_2[[1]][length(unlist(split2_2))]
    split3
}

.permutation_modified<-function(mm, Par_score, RP_score)
{
    RP_score_adj<-RP_score/max(RP_score)
    NR_sum<-sum(abs(Par_score^mm) * RP_score_adj)
    sort_Par_score<-sort(Par_score,decreasing = TRUE)
    sort_name<-names(sort_Par_score)
    sort_RP_score_adj<-RP_score_adj[sort_name]
    P_hit<-sort_RP_score_adj * (abs(sort_Par_score) ^ mm) / NR_sum
    Non_sum<-length(RP_score_adj) - sum(sort_RP_score_adj)
    P_miss<-(1 - sort_RP_score_adj) / Non_sum
    P_hit - P_miss
}
.permutation_random<-function(mm, Par_score, RP_score)
{
    RP_score_adj<-RP_score/max(RP_score)
    sort_Par_score<-sort(Par_score,decreasing=T)
    sort_name<-names(sort_Par_score)
    sort_RP_score_adj<-sample(RP_score_adj)
    NR_sum<-sum(abs(sort_Par_score^mm)* sort_RP_score_adj)
    P_hit<-sort_RP_score_adj*(abs(sort_Par_score)^mm)/NR_sum
    Non_Sum<-length(RP_score_adj)-sum(sort_RP_score_adj)
    P_miss<-(1-sort_RP_score_adj)/Non_Sum
    P_hit-P_miss
}

