all.assign.agdex <-
function(comp.vec)
{
  x <- comp.vec
	if(!is.vector(x)) x<-as.vector(x)
    	y <- sort(x)
	distinct.y <- unique(y)
	grp <- length(distinct.y)
	grp.n <- as.vector(table(y))
	L <- matrix(y[1:grp.n[1]],grp.n[1],1)
	nobs<-grp.n[1]
	for (i in 2:grp)
		for (j in 1:grp.n[i])
		{
			nobs=nobs+1
	   	L <- rbind(L,rep(distinct.y[i],dim(L)[2]))
	   	M <- choose.multinomial(nobs,as.vector(table(L[,1])))
			target.row=1
			start.row=nobs
			L1 <- L
			m <- 1
			while (m < (M-0.5))
			{
				mm <-(L1[target.row,]!=distinct.y[i])
				n.mm <- sum(mm)
				L1 <- matrix(L1[,mm],ncol=n.mm)
				temp <- L1[target.row,]
				L1[target.row,] <- L1[start.row,]
				L1[start.row,] <- temp
				L <- cbind(L,L1)
				start.row <- target.row
				target.row <- target.row+1
				m <- dim(L)[2]
			}
		}
	return(L)
}

