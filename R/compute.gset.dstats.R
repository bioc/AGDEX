compute.gset.dstats <-
function(dstats,            # vector of d-stats from one experiment
                              gset.index.list)   # list of row-indices for gene-sets

{
   n.gset <- length(gset.index.list)
   gset.dstats <- rep(NA,n.gset)
   abs.dstats <- abs(dstats)
   for (i in 1:n.gset)
   {
     gset.index <- unlist(gset.index.list[[i]])
     gset.dstats[i] <- mean(abs.dstats[gset.index])
   }
   return(gset.dstats)
}

