\name{gset.data}
\alias{gset.data}
\title{a sample gene-set data}
\description{\emph{gset.data} A gene-set data belongs to the \emph{GeneSetCollection} class of \emph{GSEABase} package. 
  
 }

\usage{
data(gset.data)
} 
\details{
This sample gene-sets data contain 10 small gene-sets which are randomly selected from the full pathway 
gene-sets downloaded from http://www.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/3.0/c5.all.v3.0.orig.gmt.
The full gene-set data contain 1454 gene-sets. For the sample gene-set data of AGDEX package, 
10 gene-sets are randomly selected, each has 20 to 30 probe-sets. Each gene-set has 20 to 30 probe-sets.  
We used \emph{getGmt} function from \emph{GSEABase} to read \emph{GMT} format into a 
\emph{GeneSetCollection} class object, then map the genes to probe-set IDs using \emph{hgu133 plus2} annotation data which contains the 
mapping from genes to probe-sets indentifiers.
}
    

\seealso{

\link[GSEABase]{GeneSetCollection}-class

\code{\link{agdex}}; \code{\link{human.data}}; \code{\link{mouse.data}}; \code{\link{map.data}}
}

\keyword{misc}

\examples{

# download the pathway gene-sets data #
\dontrun{
 gset.url <- "http://www.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/3.0/c5.all.v3.0.orig.gmt"
 gset.file.name <- unlist(strsplit(gset.url,split="/"))
 gset.file.name <- gset.file.name[length(gset.file.name)]
 gset.destination <- paste(local.data.dir,gset.file.name,sep="")
 download.file(gset.url, gset.destination)
 gset.file <- gset.destination
 gset.data <- getGmt(gset.file)        
 
# read in human U133+2 array annotation file#
 human.ann.data <- read.table("local human U133+2 array annotation data", head=T, sep="\t")
 genes.in.ann <- human.ann.data[,3]      # get the gene symbols from annotation file

# map the genes to probe-set IDs using human annotation data.########## 
 for (i in 1:length(gset.data))
  {
  genes.this.gset <- geneIds(gset.data[[i]])
  match.rows <- is.element(genes.in.ann, genes.this.gset) 
  probe.this.gset <- human.ann.data$ID[match.rows] 
  geneIds(gset.data[[i]]) <- as.character(probe.this.gset) 
  }
} 
}
