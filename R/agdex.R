agdex <- 
function(dex.setA,         # differential expression set A, list with express.set, comp.def, comp.vec, and gene-set index list
                dex.setB,         # differential expression set B,
                map.data,         # list with probe.map data frame, map.Aindex.col, map.Bindex.col to point to columns with the row indices
                min.nperms=100,   # minimum number of permutations for adaptive permutation
                max.nperms=10000) # maximum number of permutations for adaptive permutation
                         
{
   
   # Prepare data for analysis
   cat("Preparing differential expression data sets (dex.setA and dex.setB):",date(),"\n")
   dex.setA <- prep.dex.set(dex.setA)
   dex.setB <- prep.dex.set(dex.setB)
   gc()
   gc()
   do.gset  <-  !(is.null(dex.setA$gset.index.list)& is.null(dex.setB$gset.index.list))
   # extract expression matrix data from each dex.set
   exprs.A  <-  exprs(dex.setA$express.set)
   exprs.B  <-  exprs(dex.setB$express.set)
   if (is.null(rownames(exprs.A))) rownames(exprs.A) <- paste("A_probe",1:(dim(exprs.A)[1]),sep="_")
   if (is.null(rownames(exprs.B))) rownames(exprs.B) <- paste("B_probe",1:(dim(exprs.B)[1]),sep="_")
   
   A.probes <- rownames(exprs.A)
   B.probes <- rownames(exprs.B)
   
   # Prepare the map data
   cat("Preparing data that maps probe set IDs across experiments:",date(),"\n")
   map.data <- prep.map.data(map.data,A.probes,B.probes)
   
   
   #print(map.gsets)

  

   # Compute observed dstats with matrix multiplication
   a.index <- unlist(map.data$probe.map[,map.data$map.Aindex.col])
   b.index <- unlist(map.data$probe.map[,map.data$map.Bindex.col])

   cat("Computing statistics for observed data:",date(),"\n")
   t1<-proc.time()
   obs.dstatsA <- exprs.A%*%dex.setA$comp.vec
   obs.dstatsB <- exprs.B%*%dex.setB$comp.vec

   # Compute observed agdex stats
   obs.gwide.agdex.stat <- compute.agdex.stat(obs.dstatsA[a.index],
                                            obs.dstatsB[b.index])   
   
   
   obs.gset.dstatsA <- obs.gset.dstatsB <- obs.gset.agdex.stats <- map.gsets <- NULL
   if(do.gset)   
   {
    cat("Mapping gene-sets across experiments:",date(),"\n")
    map.gsets <- map.gset.lists(map.data,
                             dex.setA$gset.index.list,
                             dex.setB$gset.index.list)
   dex.setA$gset.index.list <- map.gsets$gset.listA
   dex.setB$gset.index.list <- map.gsets$gset.listB                         
   # Compute the observed gene-set AGDEX statistics
   obs.gset.agdex.stats <- compute.gset.agdex.stats(obs.dstatsA,obs.dstatsB,map.gsets$agdex.list)
   # Compute observed gset.dstats
   obs.gset.dstatsA <- compute.gset.dstats(obs.dstatsA,dex.setA$gset.index.list)
   obs.gset.dstatsB <- compute.gset.dstats(obs.dstatsB,dex.setB$gset.index.list)
   t2 <- proc.time()
   cat("Computing time for observed statistics (in seconds):","\n")
   cat(t2-t1,"\n")
   }
   

   # Now perform permutations for experiment A
   nposs.permsA <- choose(length(dex.setA$comp.vec),sum(dex.setA$comp.vec > 0)) # Determine number of possible permutations
   approx.npermsA.null <- min.nperms*(1+log(max.nperms)-log(min.nperms))      # From Pounds et al, BIBM09 conference proceedings
   
   max.npermsA <- min(max.nperms,nposs.permsA)       # Define the maximum number of permutations for experiment A
   if (nposs.permsA < approx.npermsA.null) min.npermsA <- nposs.permsA  # Might as well do exact test for everything in this case
   else min.npermsA <- min.nperms
   
   exactA <-(max.npermsA==nposs.permsA)    # Determine whether the genome-wide and gene-specific tests are exact
   if (exactA)   # Generate all possible assignments if exact test is used
   {
     asgnA <- all.assign.agdex(dex.setA$comp.vec)
     asgn.col.index <- sample(dim(asgnA)[2])  # So that the permutations are in random order for adaptive testing
     asgnA <- asgnA[,asgn.col.index]
   }

   # Initialize some variables
   ngsetA <- length(obs.gset.dstatsA)    # number of gene-sets for experiment A
   gset.keepA <- rep(T,ngsetA)           # indicate which gene-sets to keep
   gset.npermsA <- rep(0,ngsetA)         # number of permutations performed for each gene-set
   dstatsA.pval <- rep(0,length(obs.dstatsA))  # p-values for probe-set dstats
   gset.dstatsA.pval <- rep(0,ngsetA)          # p-values for gene-set dstats
   gset.agdex.pvalA <- matrix(0,ngsetA,2)           # p-values for gene-set agdex
   gwide.agdex.pvalA <- rep(0,2)
   
   cat("Permuting experiment A data",max.npermsA,"times:",date(),"\n")
   t1 <- proc.time()
   for (i in 1:max.npermsA)  # Loop over permutations
   {
      if (exactA) comp.vec <- asgnA[,i]
      else comp.vec<-sample(dex.setA$comp.vec)  # generate the permutation
      
      # Compute individual gene stats and genome-wide AGDEX stat for all permutations
      dstatsA <- exprs.A%*%comp.vec                              # Compute permutation dstat for each probe set in experiment A
      dstatsA.pval <- dstatsA.pval+(abs(dstatsA) >= abs(obs.dstatsA))   # Update p-value

      # Compute genome-wide agdex statistic for all permutations
      gwide.agdex.stat <- compute.agdex.stat(dstatsA[a.index],
                                           obs.dstatsB[b.index]) # Compute genome-wide agdex statistic
      gwide.agdex.pvalA <- gwide.agdex.pvalA+(abs(gwide.agdex.stat) >= abs(obs.gwide.agdex.stat))         # Update p-value
      #print(any(gset.keepA))
      if (any(gset.keepA))  # Compute gene-set results, if necessary
      {
         nkeep <- sum(gset.keepA)  # number of gene-sets for this permutation round
         
         gset.dstatsA <- compute.gset.dstats(dstatsA,dex.setA$gset.index.list[gset.keepA])   # Compute permutation gset.dstatsA for gset.keepA
         gset.dstatsA.pval[gset.keepA] <- gset.dstatsA.pval[gset.keepA]+
                                        (abs(gset.dstatsA) >= abs(obs.gset.dstatsA[gset.keepA])) # Update p-value for gset.keepA

        gset.agdex.stats <- compute.gset.agdex.stats(dstatsA,obs.dstatsB,map.gsets$agdex.list[gset.keepA]) # Compute permutation gset.agdex.stats for gset.keepA
        gset.agdex.pvalA[gset.keepA,] <- gset.agdex.pvalA[gset.keepA,]+
                                        (abs(gset.agdex.stats) >= abs(obs.gset.agdex.stats[gset.keepA,]))  # Update p-value
                                      
        gset.npermsA[gset.keepA] <- gset.npermsA[gset.keepA]+1   # Update permutation count

        # Update keep vector
        min.gset.pval <- apply(cbind(gset.dstatsA.pval[gset.keepA],
                                   matrix(gset.agdex.pvalA[gset.keepA,],nkeep,2)),
                             1,min)
        keep.next.round <- (min.gset.pval<min.npermsA)
        gset.keepA[gset.keepA] <- keep.next.round
      }
   }
   
   # Perform division for p-values
   dstatsA.pval <- dstatsA.pval/max.npermsA
   gwide.agdex.pvalA <- gwide.agdex.pvalA/max.npermsA
   
   gset.agdex.pvalA <- gset.agdex.pvalA/gset.npermsA
   gset.dstatsA.pval <- gset.dstatsA.pval/gset.npermsA
   t2 <- proc.time()
   cat("Computing time for permutation analysis of data set A (in seconds):","\n")
   cat(t2-t1,"\n")
   
   # gene-set results for permutation of data set A
   A.gset.res <- NULL
   if(do.gset)
   A.gset.res <- cbind.data.frame(gset.source=map.gsets$gset.source,
                                gset.name=map.gsets$gset.names,
                                A.gset.dstat=obs.gset.dstatsA,
                                A.gset.dpval=gset.dstatsA.pval,
                                A.gset.cos.stat=obs.gset.agdex.stats$cos.stat,
                                A.gset.cos.pval=gset.agdex.pvalA[,1],
                                A.gset.dop.stat=obs.gset.agdex.stats$prop.stat,
                                A.gset.dop.pval=gset.agdex.pvalA[,2],
                                A.gset.nperms=gset.npermsA)
                                
   # Now perform permutations for experiment B
   nposs.permsB <- choose(length(dex.setB$comp.vec),sum(dex.setB$comp.vec > 0)) # Determine number of possible permutations
   approx.npermsB.null <- min.nperms*(1+log(max.nperms)-log(min.nperms))      # From Pounds et al, BIBM09 conference proceedings
   
   max.npermsB <- min(max.nperms,nposs.permsB)       # Define the maximum number of permutations for experiment B
   if (nposs.permsB<approx.npermsB.null) min.npermsB <- nposs.permsB  # Might as well do exact test for everything in this case
   else min.npermsB <- min.nperms
   
   exactB <- (max.npermsB==nposs.permsB)    # Determine whether the genome-wide and gene-specific tests are exact
   if (exactB)   # Generate all possible assignments if exact test is used
   {
     asgnB <- all.assign.agdex(dex.setB$comp.vec)
     asgn.col.index <- sample(dim(asgnB)[2])  # So that the permutations are in random order for adaptive testing
     asgnB <- asgnB[,asgn.col.index]
   }

   # Initialize some variables
   ngsetB <- length(obs.gset.dstatsB)    # number of gene-sets for experiment B
   gset.keepB <- rep(T,ngsetB)           # indicate which gene-sets to keep
   gset.npermsB <- rep(0,ngsetB)         # number of permutations performed for each gene-set
   dstatsB.pval <- rep(0,length(obs.dstatsB))  # p-values for probe-set dstats
   gset.dstatsB.pval <- rep(0,ngsetB)          # p-values for gene-set dstats
   gset.agdex.pvalB <- matrix(0,ngsetB,2)           # p-values for gene-set agdex
   gwide.agdex.pvalB <- rep(0,2)
   cat("Permuting experiment B data",max.npermsB,"times:",date(),"\n")
   t1 <- proc.time()   
   for (i in 1:max.npermsB)  # Loop over permutations
   {
      if (exactB) comp.vec <- asgnB[,i]
      else comp.vec <- sample(dex.setB$comp.vec)  # generate the permutation
      
      # Compute individual gene stats and genome-wide AGDEX stat for all permutations
      dstatsB <- exprs.B%*%comp.vec                              # Compute permutation dstat for each probe set in experiment B
      dstatsB.pval <- dstatsB.pval+(abs(dstatsB) >= abs(obs.dstatsB))   # Update p-value

      # Compute genome-wide agdex statistic for all permutations
      gwide.agdex.stat <- compute.agdex.stat(dstatsB[map.data$probe.map[,map.data$map.Bindex.col]],
                                           obs.dstatsA[map.data$probe.map[,map.data$map.Aindex.col]]) # Compute genome-wide agdex statistic
      gwide.agdex.pvalB <- gwide.agdex.pvalB+(abs(gwide.agdex.stat) >= abs(obs.gwide.agdex.stat))         # Update p-value
      #print(any(gset.keepB))
      if (any(gset.keepB))  # Compute gene-set results, if necessary
      {                                               
         nkeep <- sum(gset.keepB)  # number of gene-sets for this permutation round
         
         gset.dstatsB <- compute.gset.dstats(dstatsB,dex.setB$gset.index.list[gset.keepB])   # Aompute permutation gset.dstatsB for gset.keepB
         gset.dstatsB.pval[gset.keepB] <- gset.dstatsB.pval[gset.keepB]+
                                        (abs(gset.dstatsB) >= abs(obs.gset.dstatsB[gset.keepB])) # Update p-value for gset.keepB

        gset.agdex.stats <- compute.gset.agdex.stats(obs.dstatsA,dstatsB,map.gsets$agdex.list[gset.keepB]) # Aompute permutation gset.agdex.stats for gset.keepB
        gset.agdex.pvalB[gset.keepB,] <- gset.agdex.pvalB[gset.keepB,]+
                                        (abs(gset.agdex.stats) >= abs(obs.gset.agdex.stats[gset.keepB,]))  # Update p-value
                                      
        gset.npermsB[gset.keepB] <- gset.npermsB[gset.keepB]+1   # Update permutation count

        # Update keep vector
        min.gset.pval <- apply(cbind(gset.dstatsB.pval[gset.keepB],
                                   matrix(gset.agdex.pvalB[gset.keepB,],nkeep,2)),
                             1,min)
        keep.next.round <- (min.gset.pval<min.npermsB)
        gset.keepB[gset.keepB] <- keep.next.round
      }
   }
   
   # Perform division for p-values
   dstatsB.pval <- dstatsB.pval/max.npermsB
   gwide.agdex.pvalB <- gwide.agdex.pvalB/max.npermsB
   
   gset.agdex.pvalB <- gset.agdex.pvalB/gset.npermsB
   gset.dstatsB.pval <- gset.dstatsB.pval/gset.npermsB
   
   t2<-proc.time()
   cat("Computing time for permutation analysis of data set B (in seconds):","\n")
   cat(t2-t1,"\n")
   
   # gene-set results for permutation of data set B
   B.gset.res <- gset.result <- NULL
   if(do.gset)
   {B.gset.res <- cbind.data.frame(gset.source=map.gsets$gset.source,
                                gset.name=map.gsets$gset.names,
                                B.gset.dstat=obs.gset.dstatsB,
                                B.gset.dpval=gset.dstatsB.pval,
                                B.gset.cos.stat=obs.gset.agdex.stats$cos.stat,
                                B.gset.cos.pval=gset.agdex.pvalB[,1],
                                B.gset.dop.stat=obs.gset.agdex.stats$prop.stat,
                                B.gset.dop.pval=gset.agdex.pvalB[,2],
                                B.gset.nperms=gset.npermsB)

   # Compute meta-enrichment analysis                             
   zstatA <- qnorm((A.gset.res$A.gset.nperms*A.gset.res$A.gset.dpval+0.5)/(A.gset.res$A.gset.nperms+1))*-1
   zstatB <- qnorm((B.gset.res$B.gset.nperms*B.gset.res$B.gset.dpval+0.5)/(B.gset.res$B.gset.nperms+1))*-1
   meta.zstat <- (zstatA+zstatB)/sqrt(2)
   meta.pval <- pnorm(-meta.zstat)   
          
   cat("Packaging Result Object:",date(),"\n")
   gset.result <- cbind.data.frame(A.gset.res,B.gset.res[,3:9],
                                 meta.enrich.zstat=meta.zstat,
                                 meta.enrich.pval=meta.pval)
   }
   gwide.agdex.result <- cbind.data.frame(stat.name=c("cos","dop"),
                                  stat.value=obs.gwide.agdex.stat,
                                  A.pval=gwide.agdex.pvalA,
                                  B.pval=gwide.agdex.pvalB,
                                  A.nperms=max.npermsA,
                                  A.exact=exactA,
                                  B.nperms=max.npermsB,
                                  B.exact=exactB)
                                  
   dex.resA <- cbind.data.frame(probe.id=rownames(exprs.A),
                              dstat=obs.dstatsA,
                              dpval=dstatsA.pval)
                              
   dex.resB <- cbind.data.frame(probe.id=rownames(exprs.B),
                              dstat=obs.dstatsB,
                              dpval=dstatsB.pval)
               
   # Meta-analysis for matched genes     
   one.sided.pvalA <- (gwide.agdex.result$A.nperms[1]*dex.resA$dpval+0.5)/
                    (2*gwide.agdex.result$A.nperms[1]+1)

   meta.dex.zstatA <- abs(qnorm(one.sided.pvalA))*sign(dex.resA$dstat)
   meta.dex.zstatB <- sign(dex.resB$dstat)*abs(qnorm((gwide.agdex.result$B.nperms[1]*dex.resB$dpval+0.5)/
                                          (2*gwide.agdex.result$B.nperms[1]+1)))

   meta.dex.zstat <- meta.dex.zstatA[map.data$probe.map[,map.data$map.Aindex.col]]+
                   meta.dex.zstatB[map.data$probe.map[,map.data$map.Bindex.col]]
   meta.dex.zstat <- meta.dex.zstat/sqrt(2)
   meta.dex.pval <- 2*pnorm(-abs(meta.dex.zstat))
   
   meta.dex.res <- cbind.data.frame(map.data$probe.map,
                                  dstatA=dex.resA$dstat[map.data$probe.map[,map.data$map.Aindex.col]],
                                  dpvalA=dex.resA$dpval[map.data$probe.map[,map.data$map.Aindex.col]],
                                  dstatB=dex.resB$dstat[map.data$probe.map[,map.data$map.Bindex.col]],
                                  dpvalB=dex.resB$dpval[map.data$probe.map[,map.data$map.Bindex.col]],
                                  meta.zstat=meta.dex.zstat,
                                  meta.pval=meta.dex.pval)
                                  
   cat("Done:",date(),"\n")                               
   return(list(dex.compA=dex.setA$comp.def,
               dex.compB=dex.setB$comp.def, 
               gwide.agdex.result=gwide.agdex.result,
               gset.result=gset.result,
               meta.dex.res=meta.dex.res,
               dex.resA=dex.resA,
               dex.resB=dex.resB,
               dex.asgnA=cbind.data.frame(id=sampleNames(dex.setA$express.set),
                                          grp.lbl=pData(dex.setA$express.set)[,dex.setA$comp.var]),
               dex.asgnB=cbind.data.frame(id=sampleNames(dex.setB$express.set),
                                          grp.lbl=pData(dex.setB$express.set)[,dex.setB$comp.var]),
                              gset.listA=map.gsets$gset.listA,
               gset.listB=map.gsets$gset.listB,
               gset.list.agdex=map.gsets$agdex.list))   

}

