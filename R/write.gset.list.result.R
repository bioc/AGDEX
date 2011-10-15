write.gset.list.result <-
function(agdex.res,out.file)

{
    ngset <- nrow(agdex.res$gset.result)
    if(!is.null(ngset))
    {gset.agdex <- gset.B <- gset.A <- matrix("NA", ngset, 2)
    gset.A[,1] <- gset.B[,1] <- gset.agdex[,1] <- agdex.res$gset.result$gset.name
    for (i in 1:ngset)
    { 
     A.probes <- agdex.res$dex.resA$probe.id
     B.probes <- agdex.res$dex.resB$probe.id  
     idxA <- unlist(agdex.res$gset.listA[[i]])
     idxB <- unlist(agdex.res$gset.listA[[i]])     
     gset.A[i,2] <- paste(A.probes[idxA],collapse=",")
     gset.B[i,2] <- paste(B.probes[idxB],collapse=",")
     idx1 <- unlist(agdex.res$gset.list.agdex[[i]][,1])
     idx2 <- unlist(agdex.res$gset.list.agdex[[i]][,1])
     gset.agdex[i,2] <- paste(paste(A.probes[idx1],B.probes[idx2],sep="&"),collapse=",")
    }
    colnames(gset.A) <- colnames(gset.B)<-c("gset.id","probe.ids")
    colnames(gset.agdex) <- c("gset.id", "probe.pair")
    }
    else gset.A <- gset.B <- gset.agdex <- NULL
    write("# Enrichment Gene-Set list for Comparison A:",out.file,append=T)
    write.table(gset.A,out.file,sep="\t",col.names=T,row.names=F,quote=T,append=T)
      
    write("# Enrichment Gene-Set List for Comparison B:",out.file,append=T)
    write.table(gset.B,out.file,sep="\t",col.names=T,row.names=F,quote=T,append=T)
   
    write("# AGDEX Gene-Set List:",out.file,append=T)
    write.table(gset.agdex,out.file,sep="\t",col.names=T,row.names=F,quote=T,append=T)
}

