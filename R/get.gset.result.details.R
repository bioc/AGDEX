get.gset.result.details <-
function(agdex.result,gset.ids=NULL,alpha=0.01)

{
   if (is.null(gset.ids)&(is.na(alpha)|(alpha>1)|(alpha<=0))) stop("Provide a non-empty gene-set list or set alpha between 0 and 1.")
 #  if (is.character(agdex.result))  agdex.result<-read.agdex.result(agdex.result)
    gset.name<-agdex.result$gset.result$gset.name
   if (is.null(gset.ids))
   {
      pvals <- agdex.result$gset.result[,c(4, 6, 8,13,15,18)]
      rs <- rowSums(pvals <= alpha)
      rs.keep <-rs > 0
      rs.keep[is.na(rs)] <- F
      if (!any(rs.keep)) stop(paste("No gene-sets have a p-value less than alpha =",alpha))
      gset.ids<-gset.name[rs.keep]
   }

   if (any(!is.element(gset.ids,gset.name))) warning("Some requested gene-sets do not exist in agdex.result object.")
   if (!any(is.element(gset.ids,gset.name))) stop("No requested gene-set ids exist in gene-set results of agdex.result object.")

   ngset <- length(gset.ids)
   gsetA.details <- NULL
   A.probes <- agdex.result$dex.resA$probe.id
   B.probes <- agdex.result$dex.resB$probe.id
   for (i in 1:ngset)
   {
      gset.mtch <- is.element(gset.name,gset.ids[i])
      gset.res <- agdex.result$gset.result[gset.mtch,c("gset.name","meta.enrich.zstat", "meta.enrich.pval")]
      #gset.res<-cbind.data.frame(gset.id=gset.ids[i],gset.res)
      rownames(gset.res)<- NULL
      idx.gs <- (1:length(gset.mtch))[gset.mtch]
      gene.mtch <- unlist(agdex.result$gset.listA[[idx.gs]])
      gene.res <- agdex.result$dex.resA[gene.mtch,]
      rownames(gene.res) <- NULL
      gset.detail<-cbind.data.frame(gene.res,gset.res)
      gsetA.details<-rbind.data.frame(gsetA.details,gset.detail)
   }

   gsetB.details <- NULL
   for (i in 1:ngset)
   {
      gset.mtch <- is.element(gset.name,gset.ids[i])
      gset.res <- agdex.result$gset.res[gset.mtch,c("gset.name","meta.enrich.zstat","meta.enrich.pval")]
      #gset.res<-cbind.data.frame(gset.id=gset.ids[i],gset.res)
      rownames(gset.res)<- NULL
      idx.gs <- (1:length(gset.mtch))[gset.mtch]
      gene.mtch <- unlist(agdex.result$gset.listB[[idx.gs]])
      gene.res <- agdex.result$dex.resB[gene.mtch,]
      rownames(gene.res) <- NULL
      gset.detail <- cbind.data.frame(gene.res, gset.res)
      gsetB.details <- rbind.data.frame(gsetB.details,gset.detail)
   }

   agdex.details<-NULL
   for (i in 1:ngset)
   {
      gset.mtch <- is.element(gset.name,gset.ids[i])
      gset.res <- agdex.result$gset.res[gset.mtch,c("gset.source","gset.name",
                                                 "A.gset.cos.stat", "A.gset.cos.pval",
                                                 "A.gset.dop.stat", "A.gset.dop.pval",
                                                 "B.gset.cos.stat", "B.gset.cos.pval",
                                                 "B.gset.dop.stat", "B.gset.dop.pval")]
      rownames(gset.res)<-NULL
      idx.gs <- (1:length(gset.mtch))[gset.mtch]
      gene.mtch <- is.element(agdex.result$meta.dex.res$A.index, unlist(agdex.result$gset.list.agdex[[idx.gs]][,1]))
      gene.res <- agdex.result$meta.dex.res[gene.mtch,]
      rownames(gene.res)<-NULL
      gset.detail <- cbind.data.frame(gset.res,gene.res)
      agdex.details <- rbind.data.frame(agdex.details,gset.detail)
   }


   return(list(enrichA.details=gsetA.details,
               enrichB.details=gsetB.details,
               agdex.details=agdex.details))
}

