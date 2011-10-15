\name{make.dex.set.object}
\alias{make.dex.set.object}
\title{Make a list object for function agdex()}
\description{
This function generates a list object containing four components for function agdex() 
} % end description
\usage{
make.dex.set.object(Eset.data, comp.var, comp.def, gset.collection)
}

\arguments{
  \item{Eset.data}{an \emph{ExpressionSet} object carries the gene expression data (Exprs) and Phenotype data (pData)}
  \item{comp.var}{the column name or numeric index for group labels in \emph{pData} of object \emph{Eset.data}}
  \item{comp.def}{a string definition of comparison, group labels connected by "-"}
  \item{gset.collection}{an object belongs to class \emph{GeneSetCollection}}

}

\details{
The ExpressionSet includes two components:
\emph{exprs}: a matrix of expression values 
\emph{pData}: a data frame contains the sample IDs and their assigned group labels.

\emph{gset.collection} contains a GeneSetCollection object defined in the Bioconductor package GSEABase.  
The \emph{gset.collection} object must use the same identifiers for probe-sets as that used in the \emph{exprs} component of \emph{Eset.data}.  
}

\value{
A list object containing the four components described in \emph{argument} section.
}

\author{
Stan Pounds <\email{stanley.pounds@stjude.org}; Cuilan Lani Gao <\email{cuilan.gao@stjude.org}>
}

\seealso{
ExpressionSet class: \link[Biobase]{ExpressionSet}.

GeneSetCollection class: \link[GSEABase]{GeneSetCollection}.

\code{\link{agdex}}; \code{\link{write.agdex.result}}; \code{\link{agdex.scatterplot}};
\code{\link{get.gset.result.details}}; \code{\link{read.agdex.gset.details}}; \code{\link{read.agdex.gset.details}}; 
}

\examples{
 # load data
 data(human.data)
 data(mouse.data) 
 data(gset.data)                  
 # make dex.set object for human data
 dex.set.human <- make.dex.set.object(human.data,
                                      comp.var=2,
                                      comp.def="human.tumor.typeD-other.human.tumors",
                                      gset.collection=gset.data)
 # make dex.set object for mouse data
 dex.set.mouse <- make.dex.set.object(mouse.data,
                                      comp.var=2,
                                      comp.def="mouse.tumor-mouse.control",
                                      gset.collection=NULL)
  
}




