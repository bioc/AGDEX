\name{read.agdex.gset.details}
\alias{read.agdex.gset.details}
\title{read detailed results of gene-set analyses}
\description{
A function to read output file saved by write.agdex.gset.details.
}
\usage{
read.agdex.gset.details(gset.detail.file)
}
\arguments{
  \item{gset.detail.file}{
the file name of the probe-set level gene-set result saved by function \emph{write.agdex.gset.details}
}
}

\author{
Stan Pounds<\email{stanley.pounds@stjude.org}; Cuilan Lani Gao<\email{cuilan.gao@stjude.org}>
}

\seealso{
\code{\link{agdex}};\code{\link{write.agdex.gset.details}}; \code{\link{get.gset.result.details}}
}
\examples{
#read agdex gene-set details from an output .txt file written by function "write.agdex.gset.details"
\dontrun{
 read.agdex.gset.details("gset.details.txt")
} 
}

