read.enrich.gset.list <-
function(file.name,skip,nrows)

{
   temp <- read.table(file.name,sep="\t",as.is=T,header=T,skip=skip,nrows=nrows,quote='"')
   ngset <- dim(temp)[1]
   res <- vector("list",ngset)
   names(res) <- temp[,1]
   for (i in 1:ngset) res[[i]] <- unlist(strsplit(temp[i,2],split=","))
   return(res)
}

