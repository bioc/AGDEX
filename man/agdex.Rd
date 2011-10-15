\name{agdex}
\alias{agdex}
\title{Agreement of Differential Expression Analysis}
\description{
This function performs agreement of differential expression (AGDEX) analysis across a pair of two-group experiments. 
AGDEX measures and determines the statistical 
significance of the similarity of the results from two experiments that measure differential expression across two groups.  
A metric of agreement is defined to measure the similarity and the significance is determined by permutation of group labels.
 Please see our methodology paper for details [1] (Pounds et al. 2011). 
} 
\usage{
agdex(dex.setA, dex.setB, map.data, min.nperms = 100, max.nperms = 10000)
}

\arguments{
  \item{dex.setA}{A list object with 4 components that defines a two-group comparison "A", for example "human tumor-human control". These components are \emph{express.set}, \emph{comp.def}, \emph{comp.variable} and 
  \emph{gset.collection} (optional). The \emph{express.set} component is a Bioconductor \emph{ExpressionSet} object with a matrix of expression 
  data in \emph{exprs} 
   and the phenotype data in \emph{pData}. The \emph{comp.variable} component gives the name or numeric index of the column of group label in \emph{pData} 
   of \emph{express.set} object.
  \emph{comp.def} is a string with the format "tumor-control" to define a comparison of
   expression between samples labeled as "tumor" and samples labeled as "control". The \emph{gset.collection} (optional) 
   belongs to \emph{GeneSetCollection} class. See details.}
   
   
  \item{dex.setB}{A list object that defines the other two-group comparison. It has the same structure as dex.setA.}
  
  \item{map.data}{a list object with 3 components that defines how probe-sets from \emph{dex.setA} are matched with probe-sets from \emph{dex.setB}.  
  The \emph{probe.map} component is a data.frame with each row defining how probe-sets are matched across the pair of two-group comparisons. The components
  \emph{map.Aprobe.col} and \emph{map.Bprobe.col} give the names or numeric index of the column containing probe-set identifiers in \emph{dex.setA}
  and \emph{dex.setB} respectively.}

  \item{min.nperms}{minimum number of permutations for adaptive permutation testing of gene-set level results, default is set to 100.  Adaptive permutation testing permutes
                    data until observing min.nperms statistics that exceed the observed statistic in absolute value or until max.nperms permutations are performed.  Adaptive
                    permutation testing greatly reduces computational effort for permutation analysis in many genomics applications. 
                    See [2] (Pounds et al. 2011) for more details.}
   
  \item{max.nperms}{maximum number of permutations for adaptive permutation testing of gene-set level results and fixed total number of permutations for classical permutation testing of
  probe-set level results and genome-wide agreement of differential expression, default is set to 10000.}
}

\details{
Object \emph{express.set} belongs to \emph{ExpressionSet} class. \emph{express.set} includes two components:
\emph{exprs}: a matrix of gene expression data with row of probe-sets and columns of subjects.
\emph{pData}: a data frame with each row representing a sample and two columns are sample ID and sample group label. 


\emph{gset.collection} component contains a \emph{GeneSetCollection} object defined in the Bioconductor package \emph{GSEABase}.  
The \emph{gset.collection} object must be the same identifiers for probe-sets as those used in expression matrix in \emph{express.set}.  
}

\value{
A list object with the following components
\item{dex.compA}{this string echoes the \emph{comp.def} component of \emph{dex.setA} that defines the definition for two-group comparison "A",
 for example "human tumor-human control"}

\item{dex.compB}{this string echoes the \emph{comp.def} component of \emph{dex.setB} that defines the definition for two-group comparison "B".
for example "mouse tumor-mouse control"}

\item{gwide.agdex.res}{a data.frame with the agreement statistics, p-values, and number of permutations for genome-wide agreement of differential expression analysis.}
\item{gset.res}{a data.frame with results of gene-set differential expression analysis for each comparison and gene-set agreement of differential expression analysis results.}
\item{meta.dex.res}{a data.frame with results for probe-sets matched across comparisons "A" and "B".  The data.frame includes the differential expression statistic and p-value
                    from each comparison and the meta-analysis z-statistic and p-value for differential expression.}
\item{dex.resA}{a data.frame with differential expression analysis results for individual probe-sets for two-group comparison "A". The data.frame includes the probe-set identifier,
                difference of mean log-expression statistic, and the p-value.}
\item{dex.resB}{a data.frame with the same structure as \emph{dex.resA} that gives the results for two-group comparison "B".}
\item{dex.asgnA}{a data.frame that echoes the group label assignments for comparison "A"}
\item{dex.asgnB}{a data.frame that echoes the group label assignments for comparison "B"}
\item{gset.listA}{a data.frame with gene-set lists for comparison "A".  Each row indicates an assignment of a probe-set identifier to a gene-set.}
\item{gset.listB}{a data.frame with gene-set lists for comparison "B".}
\item{gset.list.agdex}{a data.frame that assigns probe-set pairs (probe-sets from comparisons A and B that query the same gene) to gene-sets for gene-set
                       agreement of differential expression analysis.}
}

\references{
 1. S.Pounds, C.Gao, R.Johnson, K.Wright, H.Poppleton, D.Finkelstein, S.leary and R.Gilbertson (2011). 
 A procedure to statistically evaluate agreement of differential expression for cross-species genomics. 
 \emph{Bioinformatics doi: 10.1093/bioinformatics/btr362}(2011).
 
 
 2. S.Pounds, X.Cao, C.Cheng, J.Yang, D. Campana, WE.Evans, C-H.Pui, and MV. Relling(2011)  
 Integrated Analysis of Pharmacokinetic, Clinical, and SNP Microarray Data using Projection onto the Most Interesting Statistical Evidence with Adaptive Permutation Testing, 
 \emph{International Journal of Data Mining and Bioinformatics}, 5:143-157.   
}

\author{
Stan Pounds <\email{stanley.pounds@stjude.org}; Cuilan Lani Gao <\email{cuilan.gao@stjude.org}>
}

\seealso{
ExpressionSet class: \link[Biobase]{ExpressionSet}.

GeneSetCollection class: \link[GSEABase]{GeneSetCollection}.

\code{\link{human.data}}; \code{\link{mouse.data}}; \code{\link{map.data}}; \code{\link{gset.data}}
\code{\link{read.agdex.result}}; \code{\link{write.agdex.result}}; \code{\link{agdex.scatterplot}}; \code{\link{get.gset.result.details}}; 
\code{\link{write.agdex.gset.details}}; \code{\link{read.agdex.gset.details}} 
}

\examples{
 # load data
 data(human.data)
 data(mouse.data)
 data(map.data)   
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
                   
 # call agdex routine
 res <- agdex(dex.set.human,dex.set.mouse,map.data,min.nperms=5,max.nperms=10)
 
 # see visualization result of the whole genome
 agdex.scatterplot(res, gset.id=NULL)
 
 # see visualization result of a specific gene-set
 agdex.scatterplot(res, gset.id="DNA_CATABOLIC_PROCESS")
 
 # get the gene-set result of a specific gene-set
 gset.detail <- get.gset.result.details(res, gset.ids="DNA_CATABOLIC_PROCESS", alpha=0.01)
  
}




