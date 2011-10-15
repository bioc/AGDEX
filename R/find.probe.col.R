find.probe.col <-
function(ex.ds.probes,probe.map)

{
   k <- dim(probe.map)[2]
   n.mtch <- rep(0,k)
   for (i in 1:k) n.mtch[i] <- sum(is.element(probe.map[,i],ex.ds.probes))
   return(min((1:k)[n.mtch==max(n.mtch)]))
}

