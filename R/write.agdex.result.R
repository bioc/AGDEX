write.agdex.result <-
function(agdex.res,out.file)

{
   write("# AGDEX Analysis Results",out.file)
   write("# Version: SJAP01",out.file,append=T)
   write(paste("# Date Generated:",date()),out.file,append=T)
   write(paste("# Comparison A Definition:",agdex.res$dex.compA),out.file,append=T)
   write(paste("# Comparison B Definition:",agdex.res$dex.compB),out.file,append=T)

   write("# Genome-Wide AGDEX Results:",out.file,append=T)
   write.table(agdex.res$gwide.agdex.res,out.file,sep="\t",col.names=T,row.names=F,append=T)

 
   write("# Gene-Set Results:",out.file,append=T)
   write.table(agdex.res$gset.result, out.file, sep="\t",col.names=T, row.names=F, append=T)
                                        
   write("# Individual Matched-Gene Results:",out.file,append=T)
   meta.res <- agdex.res$meta.dex.res
   meta.res <- meta.res[order(meta.res$meta.pval),]
   write.table(meta.res,out.file,sep="\t",col.names=T,row.names=F,append=T)
   
   write("# Individual Gene Results for Comparison A:",out.file,append=T)
   #hdrA<-matrix(c("probe.id",names(agdex.res$dex.resA)),1,1+dim(agdex.res$dex.resA)[2])
   #write.table(hdrA,out.file,sep="\t",col.names=F,row.names=F,append=T)
   dex.resA <- agdex.res$dex.resA[order(agdex.res$dex.resA$dpval),]
   write.table(dex.resA,out.file,col.names=T,row.names=F,append=T,sep="\t")
   
   write("# Individual Gene Results for Comparison B:",out.file,append=T)
   #hdrB<-matrix(c("probe.id",names(agdex.res$dex.resB)),1,1+dim(agdex.res$dex.resB)[2])
   #write.table(hdrB,out.file,col.names=F,row.names=F,append=T,sep="\t")
   dex.resA <- agdex.res$dex.resB[order(agdex.res$dex.resA$dpval),]
   write.table(agdex.res$dex.resB,out.file,col.names=T,row.names=F,append=T,sep="\t")
   
   write("# Sample Assignments for Comparison A:",out.file,append=T)
   write.table(agdex.res$dex.asgnA,out.file,sep="\t",row.names=F,col.names=T,append=T)
   
   write("# Sample Assignments for Comparison B:",out.file,append=T)
   write.table(agdex.res$dex.asgnB,out.file,sep="\t",col.names=T,row.names=F,append=T)
   
   
   write.gset.list.result(agdex.res, out.file)
   cat("Writing result completed.", "\n")
}

