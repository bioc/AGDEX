read.gset.collection <-
function(gset.file,   # name of the gene set file
                               gset.col=3,    # col. no of the gene set name in gset data
                               probe.col=1)   # col. no of the probe IDs  in gset data
{
 gset.data <- read.table(gset.file,sep="\t",as.is=T,header=T)
 gset.name <- unique(gset.data[,gset.col])
 ngset <- length(gset.name)
 m <- nrow(gset.data)
 gc.set <- vector("list", ngset)
 
 for(i in 1:ngset)
  {
   this.gset <- gset.name[i]
   idx.probe <- (1:m)[gset.data[,gset.col]==this.gset]
   probes.this.gset <- unique(gset.data[,probe.col][idx.probe])
   gs0<-GeneSet(setName=this.gset)
   geneIds(gs0) <- probes.this.gset
   gc.set[i] <- gs0
   }
   gc.set <- GeneSetCollection(gc.set)
   return(gc.set)
}

