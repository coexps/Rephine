cnvadjust <- function(object, ...)
UseMethod("cnvadjust")
# calculate the adjust expression by removing the confounding effect of CNVs through linear regression model
cnvadjust.Rephine<-function(object, ...)
{
  exp<-object$exp
  cnv<-object$cnv
  fun<-function(x)
  {
    slope <- lm(unlist(exp[x,])~unlist(cnv[x,]))[[1]][2]
    slope
  }
  cnv_lm<-apply(as.matrix(rownames(exp)),1,fun)
  cnv_factor<-cnv * cnv_lm
  reg_factor<-exp - cnv_factor
  reg_factor
}
