\name{write.agdex.gset.details}
\alias{write.agdex.gset.details}
\title{write the output of  get.gset.result.details to a tab-delimited text file.}
\description{A function to write an output file with detailed AGDEX gene-set result.}
\usage{
write.agdex.gset.details(gset.details, out.file)
}
\arguments{
  \item{gset.details}{
result produced by \emph{get.gset.result.details}
}
  \item{out.file}{name of the output file}
}    

\value{
the path and name of output file of the result of gene-level agdex analysis
}
\author{
Stan Pounds <\email{stanley.pounds@stjude.org}; Cuilan Lani Gao <\email{cuilan.gao@stjude.org}>
}

\seealso{
\code{\link{agdex}}; \code{\link{get.gset.result.details}}; \code{\link{read.agdex.gset.details}}
}
\examples{
# load the saved result run agdex routine 
data(agdex.res)
                  
# obtain all gene-set result
gset.res.all <- get.gset.result.details(agdex.res,gset.ids = NULL, alpha=0.01)

# obtain the gene set result of memember genes
gset.res0 <- get.gset.result.details(agdex.res,gset.ids = "DNA_CATABOLIC_PROCESS", alpha=0.01)

# write the gene set details to an output file
\dontrun{ 
 write.agdex.gset.details(gset.res0, "gset.details.txt")
 }
}