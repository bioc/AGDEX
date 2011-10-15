make.dex.set.object <- function(Eset.data,           # expressionSet object
                              comp.var,              # no of column of group labels in Pdata in expressionSet
                              comp.def,              # a string definition of comparison, group labels connected by "-"
                              gset.collection=NULL)  # an object of GeneSetCollection
{
 dex.set <- list(express.set=Eset.data,
                  comp.var=comp.var,
                  comp.def=comp.def,
                  gset.collection=gset.collection)
 return(dex.set)
 }