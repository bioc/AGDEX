map.gset.lists <-
function(map.data,    # list object with probe.map data frame or matrix, map.Bindex.col, map.Aindex.col to point to the columns with row indices
                         gset.listA,  # gene set list of row indices for experiment A
                         gset.listB)  # gene set list of row indices for experiment B

{
  # Determine size of the result gene set lists
  ngsetA <- length(gset.listA)
  ngsetB <- length(gset.listB)
  totn.gset <- ngsetA+ngsetB

  # Initialize the result gene set lists
  new.listA <- vector("list",totn.gset)
  new.listB <- vector("list",totn.gset)
  agdex.list <- vector("list",totn.gset)
  
  if (is.null(names(gset.listA))&!is.null(gset.listA)) names(gset.listA) <- paste("A",1:ngsetA,sep="")
  if (is.null(names(gset.listB))&!is.null(gset.listB)) names(gset.listB) <- paste("B",1:ngsetB,sep="")
  
  gset.namesA <- gset.namesB <- NULL
  
  # Fill the gene set lists for experiment A
  if (ngsetA > 0)
  {
   for (i in 1:ngsetA)
   {
     new.listA[[i]] <- gset.listA[[i]]  # same as input
     #print(head(new.listA[[i]]))
     map.rows <- is.element(map.data$probe.map[,map.data$map.Aindex.col],gset.listA[[i]])  # find matching rows in the probe.map
     #print(table(map.rows))
     agdex.list[[i]] <- matrix(map.data$probe.map[map.rows,c(map.data$map.Aindex.col,map.data$map.Bindex.col)],ncol=2)  # agdex row index matrix
     #if (i==1) print(agdex.list[[i]])
     new.listB[[i]] <- unique(agdex.list[[i]][,2]) # row indices for experiment B determined by the map
   }
   gset.namesA <- names(gset.listA)
  }
  
  # Fill the gene set lists for experiment B
  if (ngsetB > 0)
  {
   for (i in 1:ngsetB)
   {
     new.listB[[ngsetA+i]] <- gset.listB[[i]] # same as the input
     map.rows <- is.element(map.data$probe.map[,map.data$map.Bindex.col],gset.listB[[i]])    # matching rows in the probe.map
     agdex.list[[ngsetA+i]] <- map.data$probe.map[map.rows,c(map.data$map.Aindex.col,map.data$map.Bindex.col)] # agdex row index matrix
     new.listA[[ngsetA+i]] <- unique(agdex.list[[i]][,1]) # row indices for experiment A determined by the map
   }
   gset.namesB <- names(gset.listB)
  }


  return(list(gset.listA=new.listA,gset.listB=new.listB,agdex.list=agdex.list,
              gset.names=c(gset.namesA,gset.namesB),
              gset.source=c(rep("A",ngsetA),rep("B",ngsetB))))
}

