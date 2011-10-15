gsc.to.index.list <-
function(gset.collection,probe.ids)  # Modified by stan, need only rownames from expression sets
 
 {
   ngset <- length(gset.collection)  # No. of gene sets
   nprobes <- length(probe.ids)      # No. of probe sets
   probe.index <- 1:nprobes          # vector of probe indices
   res<-vector("list",ngset)       # initialize the result
   
   for (i in 1:ngset)  # loop over gene-sets
   {
      gset.probes <- geneIds(gset.collection[[i]])  # extract the probe-set ids for gene-set i
      in.gset <- is.element(probe.ids,gset.probes)  # match agains the probe.ids vector
      res[[i]] <- unlist(probe.index[in.gset])              # assign to entry i of the new list
   }
   names(res) <- names(gset.collection) # pass along the gene-set names
   return(res)                        # return results
}

