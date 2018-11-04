Rephine <- function(drugs=matrix(0,0,0), exp=matrix(0,0,0), snv=matrix(0,0,0), cnv=matrix(0,0,0), tissue=matrix(0,0,0), cellline=NULL, genes=NULL, drug=NULL, remove.zeros=FALSE)
# Construct Raphine object, with checking
# Created 27 May 2018
{
# Check exp
  exp <- as.matrix(exp)
  nsample <- ncol(exp)
  ntags <- nrow(exp)
  if (nsample>0L && is.null(colnames(exp))) colnames(exp) <- paste0("Sample", 1L:nsample)
  if (ntags>0L && is.null(colnames(exp))) rownames(exp) <- 1L:ntags

# Make object
  x <- new("Rephine", list(drugs=as.matrix(drugs), exp=as.matrix(exp), snv=as.matrix(snv), cnv=as.matrix(cnv),tissue=as.matrix(tissue)))
  if(!is.null(genes)) {
    genes <- as.data.frame(genes, stringAsFactors=FALSE)
    if(nrow(genes) != ntags) stop("Exp and genes have different number of rows")
    row.names(genes) <- row.names(exp)
    x$genes <- genes
  }

  x

}

