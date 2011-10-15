choose.multinomial <-
function(n,k)

{
	logres <- lfactorial(n)-sum(lfactorial(k))
	res <- exp(logres)
	return(res)
}

