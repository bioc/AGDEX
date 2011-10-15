read.agdex.gset.details <-
function(gset.detail.file)

{
   gset.scan <- scan(gset.detail.file,what=character(),sep="\n")
   m <- length(gset.scan)
   char1 <- substring(gset.scan,1,1)
   comm.row <- (1:m)[char1=="#"]

   enrichA <- read.table(gset.detail.file,sep="\t",header=T,quote='"',
                       as.is=T,skip=comm.row[1],nrows=comm.row[2]-comm.row[1]-2)

   enrichB <- read.table(gset.detail.file,sep="\t",header=T,quote='"',
                       as.is=T,skip=comm.row[2],nrows=comm.row[3]-comm.row[2]-2)

   agdex.gset <- read.table(gset.detail.file,sep="\t",header=T,quote='"',
                          as.is=T,skip=comm.row[3])

   res <- list(enrichA.details=enrichA,enrichB.details=enrichB,
             agdex.details=agdex.gset)

   return(res)
}

