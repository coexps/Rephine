snvselector <- function(object, ...)
UseMethod("snvselector")

snvselector.Rephine <- function(object, p=0.1, ...)
{
  drugs <- object$drugs
  snv <- object$snv
  feature_matrix<-matrix(0,ncol=ncol(drugs),nrow=nrow(snv))
  rownames(feature_matrix)<-rownames(snv)
  colnames(feature_matrix)<-colnames(drugs)

# adaptive lasso to select drug response related mutations
  coef_m <- apply(drugs,2,function(x) .coef_m(x,snv))

# use likelihood ratio test to evaluate the significance of univariate and multivariate association of selected mutation
  feature_out<-c()
  feature_filter_out<-c()
  for(i in colnames(drugs))
  {
    lrt_res <- c()
    predictor <- rownames(coef_m)[which(coef_m[,i] != 0)]
  # likelihood ratio test univariate and leave-one-out multivairate effect.
    if(length(predictor) == 0)
    {
      lrt_res <- c()
    }
    else
    {
      picked<-t(snv)[,predictor]
      fm1 <- glm(drugs[,i] ~ picked)

      if (length(predictor) == 1)
      {
        fm2<-glm(drugs[,i]~NULL)
        loo<-lrtest(fm2,fm1)
        fm3<-glm(drugs[,i]~picked)
        single<-lrtest(fm3)
        lrt_1<-c(loo[1,2]-loo[2,3],loo[2,5],single[1,2] - single[2,2], single[2,5], single[2,3])
        lrt_1<-data.frame(t(lrt_1),row.names = predictor)
        lrt_res<-rbind(lrt_res,lrt_1)
      }
      else
      {
        for(j in 1:length(predictor)) #select one mutant in turn
        {
          fm2<-glm(drugs[,i]~picked[,-j])
          fm3<-glm(drugs[,i]~picked[,j])
          loo<-lrtest(fm2,fm1) #lrt of model of all selected mutants vs model of leave one mutant out
          single<-lrtest(fm3) #lrt model of univariate
          lrt_1<-c(loo[1,2]-loo[2,2],loo[2,5],single[1,2]-single[2,2],single[2,5],single[2,3])
          lrt_res<-rbind(lrt_res,lrt_1)
        }
        rownames(lrt_res)<-predictor
      }
      feature_out<-c(feature_out, list(lrt_res))
      keep<-apply(lrt_res,1,function(x) .fun_filter(x,p))
      feature_filter_out<-c(feature_filter_out,list(names(which(keep))))
      feature_matrix[names(which(keep)),i] <- 1
    }
    # filter drug response related mutations based on the P-value of the likelihood ratio test
  }
  #feature_filter_out_all<-as.vector(unlist(feature_filter_out))
  feature_matrix<-feature_matrix[rowSums(feature_matrix) > 0,]
  feature_matrix
}

# Apply function for adaptive lasso to select drug response related mutations
.coef_m <- function(x,snv)
{
  x_var <- snv
  y_var <- x
  betals <- ginv(x_var%*% t(x_var))%*%x_var%*%y_var
  w_alasso <- 1/abs(betals)
  res2 <- cv.glmnet(t(x_var),y_var, penalty.factor = w_alasso, nfolds =  nrow(x_var),grouped = FALSE)
  coef1 <- as.matrix(coef(res2, s = "lambda.min"))[-1,]
  coef1
}

.fun_filter <- function(x,p)
{
  if(x[1] < 0 & x[2] < p & x[3] > 0 & x[4] < p)
  {
    TRUE
  }
  else
  {
    FALSE
  }
}
