\name{map.data}
\alias{map.data}
\title{Probe-set Mapping Data}
\description{A mapping between Human probe-sets and Mouse probe-sets}

\usage{
data(map.data)
}

\details{\emph{map.data} is a list object containing 3 components, \emph{probe.map}, \emph{map.Aprobe.col} and \emph{map.Bprobe.col}. 
The component \emph{probe.map} is a data frame with 2 columns of probe-sets identifiers from human and mouse array respectively. 
This sample mapping data \emph{probe.map} is a subset of the full mapping data set across human genome U133 Plus 2.0 Array and 
 mouse Expression 430 Array available at www.affymetrix.com.
  They provide mutiple match mode data sets, such as \emph{Good Match}, \emph{Complex Match} \emph{Best match} etc.
  We downloaded the \emph{Best Match} data set as our mapping data. 
  The probe-sets in \emph{probe.map} are selected such that they are contained 
  in both expression matrix of \emph{human.data} and \emph{mouse.data}.
Users can choose mapping data according to the species and platforms of their gene expression profiles either by downloading from
 www.affymetrix.com or from other sources. The array platforms of mapping data must match that of gene expression profile of 
 each species.    
    
 \emph{probe.map} component of the sample data \emph{map.data} contains 490 rows of ortholog-matched probe-sets across human array and mouse array.
 \emph{map.Aprobe.col} and \emph{map.Bprobe.col} specify column number or name containing the probe-sets from study A and study B respectively.
}

\value{
\item{probe.map}{A data frame with 490 rows and 2 columns of probe-sets identifiers of human and mouse}
\item{map.Aprobe.col}{Column number or name in data frame \emph{probe.map} containing probe-sets IDs from study A}
\item{map.Bprobel.col}{Column number or name in data frame \emph{probe.map} containing probe-sets IDs from study B} 

}

\seealso{
\code{\link{agdex}}; \code{\link{human.data}}; \code{\link{mouse.data}}; \code{\link{gset.data}}
}

\examples{

# download the "best match" mapping data across human array and mouse array # 
\dontrun{
 map.url <- "http://www.affymetrix.com/Auth/analysis/downloads/na31/ivt/Mouse430_2.na31.ortholog.csv.zip"
 map.file.name <- unlist(strsplit(map.url,split="/"))
 map.file.name <- map.file.name[length(map.file.name)]
 map.destination <- paste(local.data.dir,map.file.name,sep="")
 download.file(map.url,map.destination)  # Affy website may need users to register first
 unzip(map.destination)
 affy.ortho.file <- substring(map.destination,1,nchar(map.destination)-4)
 

# read in the mapping data  #
 ortho.data <- read.csv(affy.ortho.file,quote='"',as.is=T)

 keep the probe-sets identifers of human array "HG-U133_Plus_2" only 
 keep.rows <- is.element(ortho.data$Ortholog.Array,"HG-U133_Plus_2")
 ortho.data <- ortho.data[keep.rows,]

 ortho.data <- ortho.data[,c(1,3)]             # keep the columns containg probe-sets only
 ortho.data$Ortholog.Probe.Set <- tolower(ortho.data$Ortholog.Probe.Set)
 
# prepare the list obejct of map data for calling AGDEX routine #
 map.data <- list(probe.map=ortho.data,
                  map.Aprobe.col=2,   # the column index containing human probe-sets IDs in data frame probe.map 
                  map.Bprobe.col=1)   # the column index containing mouse probe-sets IDs in data frame probe.map
}
}


