\name{mouse.data}
\alias{mouse.data}
\title{Mouse Data}

\description{\emph{mouse.data} is an \emph{ExpressionSet} object where \emph{exprs} slot 
 carries the mouse gene expression data and \emph{pData} carries the phenotype data.}

\usage{
data(mouse.data)
}     
  
\details{
This \emph{mouse.data} is a subset taken from our published mouse data (Johnson et al. 2010).
   The original full mouse data was profiled by affymetrix 430 v2(mRNA). 
   It contains 45037 probe-sets for 13 mice brain tumors and 179 normal mice stem cells. 
 
   We used the best-match data (available from www.affymetrix.com) as the mapping data of
   ortholog-matched probe-sets across human gene expression data 
  and mouse gene expression data. \emph{expr} slot of \emph{mouse.data} 
  is a matrix of subset of the full mouse gene expression data. Those selected mouse probe-sets are ortholog-matched with 
  human probe-sets in \emph{human.data}. pData slot of the \emph{mouse.data} is a data frame with row representing mouse samples and two columns
   indicating sample IDs and sample group labels for each sample (either "mouse.tumor" or "mouse.control").

}

\value{
\item{expr(mouse.data)}{A matrix with 264 rows and 192 columns with rows representing probe-sets and cloumns of mouse
sample IDs. Each row name of the matrix is a probe-set Identifier.}

\item{pData(mouse.data)}{A data frame with 192 rows and 2 columns. Each row represents one mouse sample. Column \emph{id} is the mouse sample ID and 
\emph{grp} is the assigned sample group label.}
}

\seealso{
\link[Biobase]{ExpressionSet}-class


\code{\link{agdex}}; \code{\link{human.data}}; \code{\link{map.data}}; \code{\link{gset.data}}
}

\references{
	R. Johnson et al.(2010) Cross-species genomics matches driver mutations and cell compartments to model ependymoma. \emph{Nature}, 466: 632-6.
}

\examples{
data(mouse.data)
mouse.express.set <- exprs(mouse.data)
mouse.pheno.data <- pData(mouse.data)
}


