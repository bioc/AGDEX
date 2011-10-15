\name{human.data}
\alias{human.data}
\title{Sample ExpressionSet object of human data}
\description{An \emph{ExpressionSet} object of human data.}
   
  
\usage{
data(human.data)
}   
  
\details{
   \emph{human.data} is an \emph{ExpressionSet} object where \emph{exprs} slot 
 carries the human gene expression data and the \emph{pData} contains the phenotype data. 
 This sample data \emph{human.data} is a subset taken from our published study of human brain tumor ependymoma (Johnson et al. 2010).
   The original full human expression data contains 54,613 probe-sets for 83 human enpendymoma
   tumors. The gene expression is profiled with Affymetrix U133+2(mRNA) array and the expression data were normalized with MAS 5.0 algorithm.  
   
    
   The \emph{expr} of \emph{human.data} is a gene expression matrix with rows of probe-sets and columns representing ependymoma tumors which 
   are classified as belonging to the novel subgroup D or others. Probe-sets in the gene expression matrix are randomly selected 
   from the full human gene profile such that the selected probe-sets belong to 
    the gene-sets in \emph{gset.data}. pData slot of \emph{human.data} is a data frame 
   with two columns indicating sample ID and sample group label for each sample(either "human.tumor.typeD" or "other.human.tumors").  
}

\value{
\item{expr(human.data)}{A matrix with 246 rows and 83 columns with rows representing probe-sets and cloumns of human sample IDs. 
}

\item{pData(human.data)}{A data frame with 83 rows and 2 columns. Each row represents one human sample. Column \emph{id} is the human sample ID and 
\emph{grp} is the assigned sample group label.} 
}

\seealso{
\link[Biobase]{ExpressionSet}-class


\code{\link{agdex}}; \code{\link{mouse.data}}; \code{\link{map.data}}; \code{\link{gset.data}}
}

\references{
	R. Johnson et al.(2010) Cross-species genomics matches driver mutations and cell compartments to model ependymoma. \emph{Nature}, 466: 632-6.
}

\examples{
data(human.data)
human.express.set <- exprs(human.data)
human.pheno.data <- pData(human.data)

}
