compute.gset.agdex.stats <-
function(dstatsA,dstatsB,agdex.gset.index.list)

{
   n.gset <- length(agdex.gset.index.list)
   cos.stat <-prop.stat <- rep(NA,n.gset)
   for (i in 1:n.gset)
   {
      indA <- unlist(agdex.gset.index.list[[i]][,1])
      indB <- unlist(agdex.gset.index.list[[i]][,2])
      gset.res <- compute.agdex.stat(dstatsA[indA],dstatsB[indB])
      cos.stat[i] <- gset.res["cos.stat"]
      prop.stat[i] <- gset.res["prop.stat"]
   }
   final.result <- cbind.data.frame(cos.stat=cos.stat,prop.stat=prop.stat)
   return(final.result)
}

