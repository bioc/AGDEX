prep.dex.set <-
function(dex.set)

{
   # Just check some basic requirements
   if (is.null(dex.set$comp.def)) stop("dex.set does not define the differential expression comparison.")
   if (regexpr("-",dex.set$comp.def)<0) stop("dex.set does not include a subtraction sign.")

   grps <- unlist(strsplit(dex.set$comp.def,split="-"))
   if (length(grps)> 2) stop("dex.set$comp.def defines more than 2 groups.")
   in.data <- is.element(grps,dex.set$express.set[[dex.set$comp.var]])
   if (any(!in.data)) stop("dex.set$comp.def defines a group not in the expression set.")
   
   if (any(!is.finite(exprs(dex.set$express.set)))) stop("Expression data contains some non-numeric or infinite values.")

   # Subset the dex.set on the subjects to be included in the analysis
   if (length(grps)==1)   keep <- rep(T,dim(pData(dex.set$express.set))[1])
   else  keep <- is.element(dex.set$express.set[[dex.set$comp.var]],grps)
   dex.set$express.set <- dex.set$express.set[,keep]

   # Define the comparison vector to compute the dstat via matrix multiplication
   if (length(grps)==1)  comp.vec<-is.element(dex.set$express.set[[dex.set$comp.var]],grps[1])-
                                   !is.element(dex.set$express.set[[dex.set$comp.var]],grps[1])
   else comp.vec<-is.element(dex.set$express.set[[dex.set$comp.var]],grps[1])-
                  is.element(dex.set$express.set[[dex.set$comp.var]],grps[2])

   comp.vec[comp.vec < 0] <- comp.vec[comp.vec < 0]/sum(comp.vec < 0)    
   comp.vec[comp.vec > 0] <- comp.vec[comp.vec > 0]/sum(comp.vec > 0)

   dex.set$comp.vec <- comp.vec

   # generate gset.index.list from gset.collection
   probe.ids <- rownames(exprs(dex.set$express.set))
   if (is.null(dex.set$gset.index.list)) 
   {
      if (is.null(dex.set$gset.collection)) dex.set$gset.index.list <- NULL
      else dex.set$gset.index.list <- gsc.to.index.list(dex.set$gset.collection,probe.ids)
   }
   
   return(dex.set)
}

