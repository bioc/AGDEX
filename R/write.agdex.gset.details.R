write.agdex.gset.details <-
function(gset.details,out.file)

{
    write("# Detailed Gene-Set Enrichment Results for Experiment A:",out.file,append=F)
    write.table(gset.details$enrichA.details,out.file,sep="\t",col.names=T,row.names=F,quote=T,append=T)
    write("# Detailed Gene-Set Enrichment Results for Experiment B:",out.file,append=T)
    write.table(gset.details$enrichB.details,out.file,sep="\t",col.names=T,row.names=F,quote=T,append=T)
    write("# Detailed AGDEX Results for Paired Probe Sets:",out.file,append=T)
    write.table(gset.details$agdex.details,out.file,sep="\t",col.names=T,row.names=F,quote=T,append=T)
}

