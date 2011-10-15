compute.agdex.stat <-
function(dstatA,dstatB)

{
   x <- dstatA
   y <- dstatB
   cos.stat <- mean(x*y)/sqrt(mean(x^2)*mean(y^2))
   prop.stat <- mean(sign(x)*sign(y))
   return(c(cos.stat=cos.stat,prop.stat=prop.stat))
}

