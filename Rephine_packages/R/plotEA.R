plotEA<-function(chip,drug,col,partial,...)
{
  sort_inter.name_chip <- intersect(rownames(partial),chip[,2])
  order_merge <- order(chip[,10],decreasing=T)
  order_genes <- chip[order_merge,2]
  uniq_gene_order <- which(!duplicated(order_genes))
  uniq_value <- chip[order_merge[uniq_gene_order],col]
  uniq_value
  names(uniq_value) <- order_genes[uniq_gene_order]
  uniq_value <- uniq_value[sort_inter.name_chip]
  perm <- .permutation_modified(1,partial[sort_inter.name_chip,drug],uniq_value)
  x <- cumsum(perm)
  plot(x,...)
  rand_max <- c()
  rand_min <- c()
  for (i in 1:1000)
  {
    cum_rand <- .permutation_random(1,partial[sort_inter.name_chip,drug],uniq_value)
    rand_max <- c(rand_max,max(cum_rand))
    rand_min <- c(rand_min,min(cum_rand))
  }
  for (i in 1:6)
  {
    sum_max <- summary(rand_max)[i]
    sum_min <- summary(rand_min)[i]
    col_h <- heat.colors((9))
    abline(a=as.matrix(sum_max),b=0,col=col_h[i],lty=1,lwd=2)
    abline(a=as.matrix(sum_min),b=0,col=col_h[i],lty=1,lwd=2)
    abline(a=0,b=0)
  }
}
