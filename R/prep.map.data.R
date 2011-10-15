prep.map.data <-
function(map.data,A.probes,B.probes)

{
   nA <- length(A.probes)
   A.index <- cbind.data.frame(A.probe=A.probes,A.index=1:nA)
   
   nB <- length(B.probes)
   B.index <- cbind.data.frame(B.probe=B.probes,B.index=1:nB)
   
   probe.map1 <- map.data$probe.map

   if (is.null(map.data$map.Aprobe.col)) map.data$map.Aprobe.col <- find.probe.col(A.probes,probe.map1)
   if (is.null(map.data$map.Bprobe.col)) map.data$map.Bprobe.col <- find.probe.col(B.probes,probe.map1)
   
   names(probe.map1)[map.data$map.Aprobe.col] <- "probeA"            # lani's input for the purpose of  finding 
   names(probe.map1)[map.data$map.Bprobe.col] <- "probeB"            #  probes of each study for writing functions
   k1 <- dim(probe.map1)[2]
   probe.map2 <- merge(probe.map1,A.index,by.x="probeA", by.y="A.probe")
   map.data$map.Aindex.col <- k1+1
   
   k2 <- dim(probe.map2)[2]
   probe.map3 <- merge(probe.map2,B.index,by.x="probeB", by.y="B.probe")
   map.data$map.Bindex.col <- k2+1
   
   map.data$probe.map <- probe.map3
   
   return(map.data)
}

