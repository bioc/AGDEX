\name{get.gset.result.details}
\alias{get.gset.result.details}
\title{Extract gene-level details from gene-set results}
\description{
A function to extract the probe-sets level details for gene-set results, i.e. obtain the differential expression statistics
for probe-sets assigned to the gene-sets.  This allows users to explore which probe-sets results contribute the most to gene-set differential
expression analysis statistics and gene-set agreement of differential expression statistics.
}
\usage{
get.gset.result.details(agdex.result, gset.ids = NULL, alpha=0.01)
}

\arguments{
  \item{agdex.result}{
agdex result returned by function \emph{agdex}
}
  \item{gset.ids}{
a vector of gene-set IDs. If NULL, the result will return gene level details for all significant gene-sets at a chosen 
significant level alpha. 
}
  \item{alpha}{significance level of gene-set, default set to 0.01}
}

\value{
This function returns a list of three components.
  \item{gsetA.details }{Gene-set details result for experiment A, including differential expression statistic 
  and p-value for each probe-set in each gene-set. Each row represents a probe-set from A. The columns give gene-set name, enrichment statistic
  and its corresponding p-values, differential expression statistics and p-values}
  \item{gsetB.details }{similar result of gene-set details for experiment B. Rows and columns give the similar information to \emph{gsetA.details}.}
  \item{agdex.details}{A data frame. Each row gives results for one probe-set pair. The columns 
 give the gene-set names, cosine statistic and difference of proportions statistic and p-values, meta statistic and its p-value.}
}
\author{
Stan Pounds <\email{stanley.pounds@stjude.org}; Cuilan Lani Gao <\email{cuilan.gao@stjude.org}>
}
\seealso{
\code{\link{write.agdex.gset.details}}; \code{\link{read.agdex.gset.details}}
}

\examples{
# Load saved result run by agdex routine               
data(agdex.res)

# obtain gene-set result
gset.res.all <- get.gset.result.details(agdex.res,gset.ids = NULL, alpha=0.01)

# obtain the detailed gene set for specified gene-sets 
gset.res0 <- get.gset.result.details(agdex.res, gset.ids=c("DNA_CATABOLIC_PROCESS","GOLGI_STACK"), alpha=0.01)
}
