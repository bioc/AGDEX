\name{read.agdex.result}
\alias{read.agdex.result}
\title{Read the output file of agdex analysis}
\description{
A function to read the output file produced by the function \emph{write.agdex.result}
}
\usage{
read.agdex.result(res.file)
}

\arguments{
  \item{res.file}{
the name of output file of AGDEX analysis
}
}

\value{
Returns a list of objects of result.
\item{dex.compA}{comparison definition for "A"}
\item{dex.compB }{comparison definition for "B"}
\item{gwide.agdex.res}{genome-wide AGDEX result}
\item{gset.res}{gene-set Results}
\item{meta.dex.res}{Individual Matched-Gene results}
\item{dex.resA}{Individual Gene Results for comparison "A"}
\item{dex.resB}{Individual Gene Results from comparison "B"}
\item{dex.asgnA}{sample assignments for comparison "A"}
\item{dex.asgnB}{aample assignments for comparison "B"}
\item{gset.listA}{gene-set lists for comparison "A"}
\item{gset.listB}{gene-set lists for comparison "B"}
\item{gset.list.agdex}{gene-set lists for agdex analysis}
}
\author{
Stan Pounds<\email{stanley.pounds@stjude.org}; Cuilan Lani Gao<\email{cuilan.gao@stjude.org}>
}
\seealso{
\code{\link{agdex}};\code{\link{write.agdex.result}}
}
\examples{
#read agdex file from output file
\dontrun{
 read.agdex.result("agdex.result.txt")
 }
}