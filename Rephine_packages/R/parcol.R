parcol <- function(object, ...)
UseMethod("parcol")
require("ppcor")
require("lmtest")
require("MASS")
# calculate the partial correlation between adjusted expression and drug response by accounting for the confoundings.
parcol.Rephine <- function(object, reg_factor=matrix(0,0,0), feature_matrix=matrix(0,0,0),tcga=NULL, ...)
{
  drugs<-object$drugs
  #exp<-object$exp
  reg_factor<-reg_factor
  snv<-object$snv
  colnames(snv)->snv_cell
  colnames(reg_factor)->reg_cell
  rownames(drugs)->drug_cell
  tissue=NULL
  try(tissue<-object$tissue, silent = T)
  #print (length(tissue))
  if(length(tissue)==0){
     cell_all<-intersect(drug_cell,intersect(snv_cell, reg_cell))
     tissue2<-NULL
  }
    else
  {
      #print (1)
    cell_all<-intersect(intersect(drug_cell,rownames(as.matrix(tissue))),intersect(snv_cell, reg_cell))
    tissue2<-as.matrix(tissue)[cell_all,]
  }
  #print (tissue2)

  stage_all_partial<-c()
  for(i in colnames(drugs))
  {
      #print (cell_all)
    snv_tmp<-t(snv)[,names(which(feature_matrix[,i]!=0))]
    confounder_matrix=cbind(as.matrix(snv_tmp)[cell_all,],tissue2)
    # print (confounder_matrix)
    if (dim(confounder_matrix)[2]>0)
        #print (confounder_matrix)
        {
            #print (confounder_matrix)
            partial<-apply(reg_factor,1,function(x) unlist(.partest(drugs[,i],x[cell_all],confounder_matrix)))
            stage_all_partial<-cbind(stage_all_partial,unlist(partial))
            #print (partial)
        }
        else
        {
            partial<-apply(reg_factor,1,function(x) cor.test(drugs[,i],x[cell_all],use="pairwise.complete.obs")$estimate)
            stage_all_partial<-cbind(stage_all_partial,unlist(partial))

        }

  }
  colnames(stage_all_partial)<-colnames(drugs)
  if (length(tcga)==0){stage_all_partial}
  else
  {
      gsub(" ","",rownames(diff))->diff2_name
      rownames(diff)<-diff2_name
      apply(diff,1,function(x) sort(abs(x[-15]))[5]<=5e-3)->tttt # if one gene is differential less than 1/3 cancer types, set weight to 0 only a few genes were filtered in this step but it is quite useful to remove the random bias
      as.numeric(tttt)->tttt_num
      names(tttt_num)<-names(tttt)
      inter_names_tcga<-intersect(names(tttt),rownames(stage_all_partial))
      print(tttt_num[inter_names_tcga])
      apply((stage_all_partial)[inter_names_tcga,],2,function(x) x* tttt_num[inter_names_tcga])->stage_all_partial
  }
  stage_all_partial
}

.partest<-function(drugs,x,confounder_matrix)
{
  cbind(drugs,x,confounder_matrix)->par_matrix
  na.omit(par_matrix)->par_matrix
  pcor.test(par_matrix[,1],par_matrix[,2],par_matrix[,3:dim(par_matrix)[2]])[1]
}




