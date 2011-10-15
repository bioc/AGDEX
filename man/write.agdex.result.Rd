\name{write.agdex.result}
\alias{write.agdex.result}
\title{
Write the AGDEX results to output file
}
\description{
 A function to write the results of an AGDEX analysis to a tab-delimited text 
 output file that can be viewed in Excel or re-imported with the function \emph{read.agdex.result}. 
}
\usage{
write.agdex.result(agdex.res, out.file)
}
\arguments{
  \item{agdex.res}{
result object produced by the \emph{agdex} function
}
  \item{out.file}{
name of the output file
}
}
\author{
Stan Pounds<\email{stanley.pounds@stjude.org}; Cuilan Lani Gao<\email{cuilan.gao@stjude.org}>
}
\seealso{
\code{\link{read.agdex.result}} \code{\link{agdex}}
}
\examples{
data(agdex.res)    

\dontrun{
 #set the wording dictionary
  setwd("localWorking dictionary")
                        
 #write the agdex result to an out file
  \dontrun{
  write.agdex.result(agdex.res, "agdex.result.txt") 
  }
  
 #read the result file stored on dist back into R
  agdex.res2 <- read.agdex.result("agdex.result.txt")
 }
}

