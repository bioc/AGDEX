\name{agdex.scatterplot}
\alias{agdex.scatterplot}
\title{scatter plot of AGDEX result}
\description{
A function to visualize the results of AGDEX analysis for the entire genome or specific gene-sets in a scatterplot.}
\usage{
agdex.scatterplot(agdex.res, gset.id = NULL)
}

\arguments{
  \item{agdex.res}{
result of an AGDEX analysis, the returned result of the function \emph{agdex}.
}
  \item{gset.id}{
a specified gene-set identifier. Default is set to NULL, which produces a visualization result of the entire genome.
}
}
\value{
Returns either a scatter plot of pairs of difference in average log-expression values for genome-wide anlysis or a specified gene-set. 
}
\author{
Stan Pounds <\email{stanley.pounds@stjude.org}; Cuilan Lani Gao <\email{cuilan.gao@stjude.org}>
}
\seealso{
\code{\link{agdex}}; \code{\link{get.gset.result.details}}
}
\examples{
data(agdex.res)
# see visualization result of the whole genome
agdex.scatterplot(agdex.res, gset.id=NULL)
# scatterplot for a specified gene-set
agdex.scatterplot(agdex.res, gset.id="DNA_CATABOLIC_PROCESS") 
}