read.agdex.gset.list <-
function(file.name,skip)

{
   temp <- read.table(file.name,sep="\t",as.is=T,header=T,skip=skip, quote='"')
   ngset <- dim(temp)[1]
   res <- vector("list",ngset)
   names(res) <- temp[,1]
   for (i in 1:ngset)
   {
      str.vec <- unlist(strsplit(temp[i,2],split=","))
      npairs <- length(str.vec)
      res[[i]] <- matrix("NA",npairs,2)
      spc.mtch <- regexpr("&",str.vec)
      res[[i]][,1] <- substring(str.vec,1,spc.mtch-1)
      res[[i]][,2] <- substring(str.vec,spc.mtch+1)
   }
   return(res)
}

