read.agdex.result <-
function(res.file)

{
   # Read first character from each row
   scn <- scan(res.file,sep="\n",what=character())
   char1 <- substring(scn,1,1)
   
   # Determine which rows begin with "#"
   n<-length(char1)
   comm.row <- (1:n)[char1=="#"]
   
   # Extract comparison definitions from the file
   dex.compA <- substring(scn[4],regexpr(":",scn[4])+2)
   dex.compB <- substring(scn[5],regexpr(":",scn[5])+2)
   
   # Extract genome-wide AGDEX result
   gwide.agdex.res <- read.table(res.file,sep="\t",as.is=T,header=T,skip=6,
                               nrows=comm.row[7]-comm.row[6]-2)
   
   # Extract Gene-Set Results
   if(comm.row[8]-comm.row[7]-2 > 0)
   gset.res<-read.table(res.file,sep="\t",as.is=T,header=T,skip=comm.row[7],
                        nrows=comm.row[8]-comm.row[7]-2,quote='"', comment.char="")
   else gset.res<-NULL                     
   # Extract Individual Matched-Gene Results
   meta.dex.res <- read.table(res.file,sep="\t",as.is=T,header=T,skip=comm.row[8],
                            nrows=comm.row[9]-comm.row[8]-2,quote='"', comment.char="")
                            
   # Extract Individual Gene Results for Comparison A
   dex.resA <- read.table(res.file,sep="\t",as.is=T,header=T,skip=comm.row[9],
                        nrows=comm.row[10]-comm.row[9]-2)
                        
   # Extract Individual Gene Results for Comparison B
   dex.resB <- read.table(res.file,sep="\t",as.is=T,header=T,skip=comm.row[10],
                        nrows=comm.row[11]-comm.row[10]-2)
                        
   # Extract Sample Assignments for Comparison A
   dex.asgnA <- read.table(res.file,sep="\t",as.is=T,header=T,skip=comm.row[11],
                         nrows=comm.row[12]-comm.row[11]-2)
                         
   # Extract Sample Assignments for Comparison B
   dex.asgnB <- read.table(res.file,sep="\t",as.is=T,header=T,skip=comm.row[12],
                         nrows=comm.row[13]-comm.row[12]-2)
   if (comm.row[15]-comm.row[13] < 5) gset.listA <- gset.listB<-gset.list.agdex<-NULL                   
   else 
   {gset.listA <- read.enrich.gset.list(res.file,skip=comm.row[13],
                                    nrows=comm.row[14]-comm.row[13]-2)
                                    
   gset.listB <- read.enrich.gset.list(res.file,skip=comm.row[14],
                                    nrows=comm.row[15]-comm.row[14]-2)
                                    
   gset.list.agdex <- read.agdex.gset.list(res.file,skip=comm.row[15])
   }
   
   
   res <- list(dex.compA=dex.compA,
             dex.compB=dex.compB,
             gwide.agdex.res=gwide.agdex.res,
             gset.res=gset.res,
             meta.dex.res=meta.dex.res,
             dex.resA=dex.resA,
             dex.resB=dex.resB,
             dex.asgnA=dex.asgnA,
             dex.asgnB=dex.asgnB,
             gset.listA=gset.listA,
             gset.listB=gset.listB,
             gset.list.agdex=gset.list.agdex)
             
   return(res)
}

