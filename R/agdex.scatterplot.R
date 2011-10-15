agdex.scatterplot <-
function (agdex.res, gset.id = NULL)
{
    x.lbl <- agdex.res$dex.compA
    y.lbl <- agdex.res$dex.compB
    if (is.null(gset.id)) {
        x <- agdex.res$meta.dex.res$dstatA
        y <- agdex.res$meta.dex.res$dstatB
        cos.stat <- agdex.res$gwide.agdex.res$stat.value[agdex.res$gwide.agdex.res$stat.name ==
            "cos"]
        cos.pvalA <- agdex.res$gwide.agdex.res$A.pval[agdex.res$gwide.agdex.res$stat.name ==
            "cos"]
        cos.pvalB <- agdex.res$gwide.agdex.res$B.pval[agdex.res$gwide.agdex.res$stat.name ==
            "cos"]
        dop.stat <- agdex.res$gwide.agdex.res$stat.value[agdex.res$gwide.agdex.res$stat.name ==
            "dop"]
        dop.pvalA <- agdex.res$gwide.agdex.res$A.pval[agdex.res$gwide.agdex.res$stat.name ==
            "dop"]
        dop.pvalB <- agdex.res$gwide.agdex.res$B.pval[agdex.res$gwide.agdex.res$stat.name ==
            "dop"]
      agr <- (sign(x) == sign(y))
    sub.lbl <- paste("cos=", round(cos.stat, 2), 
                     ", cos.pvalA=",round(cos.pvalA, 4), 
                     ", cos.pvalB=", round(cos.pvalB, 4),
                     "; dop=", round(dop.stat, 2), 
                     ", dop.pvalA=", round(dop.pvalA,4), 
                     ", dop.pvalB=", round(dop.pvalB, 4))
    plot(x, y, xlab = x.lbl, ylab = y.lbl, pch = 19, col = "gray",
        sub = sub.lbl, main = gset.id)
    points(x[agr], y[agr], col = "red", pch = 19)
    lines(range(x), c(0, 0), lwd = 5, col = "black")
    lines(c(0, 0), range(y), lwd = 5, col = "black")
    pca.res<-prcomp(cbind(x,y),center=F)
    abline(a=0,b=pca.res$rotation[2,1]/pca.res$rotation[1,1],lwd=5,col="blue")
    
    }
    else {
        if (length(gset.id) != 1)
            stop("Please provide only one gene-set identifier.")
        ugset <- agdex.res$gset.result$gset.name
        ngset <- length(ugset)
        gset.mtch <- is.element(ugset, gset.id)
        if (sum(gset.mtch) < 1)
            stop(paste(gset.id, "not found in gene-set definitions of agdex.result object."))
        if (sum(gset.mtch) > 1)
            stop(paste(gset.id, "matches more than one gene-set identifier in gene-set definitions of agdex.result object."))
        gset.num <- (1:ngset)[gset.mtch]
        A.mtch<-unlist(agdex.res$gset.list.agdex[[gset.num]][,1])
     
        mtch.row<-is.element(agdex.res$meta.dex.res$A.index, A.mtch)
        x <- agdex.res$meta.dex.res$dstatA[mtch.row]
        y <- agdex.res$meta.dex.res$dstatB[mtch.row]
        gset.row<-(1:ngset)[agdex.res$gset.res$gset.name==gset.id]
        A.gset.cos.stat <- agdex.res$gset.res[gset.row, "A.gset.cos.stat"]
        B.gset.cos.stat<-agdex.res$gset.res[gset.row, "B.gset.cos.stat"]
        cos.pvalA <- agdex.res$gset.res[gset.row, "A.gset.cos.pval"]
        cos.pvalB <- agdex.res$gset.res[gset.row, "B.gset.cos.pval"]
        A.gset.dop.stat <- agdex.res$gset.res[gset.row, "A.gset.dop.stat"]
        B.gset.dop.stat <- agdex.res$gset.res[gset.row, "B.gset.dop.stat"]
        dop.pvalA <- agdex.res$gset.res[gset.row, "A.gset.dop.pval"]
        dop.pvalB <- agdex.res$gset.res[gset.row, "B.gset.dop.pval"]
     agr <- (sign(x) == sign(y))
    sub.lbl <- paste("cos=", round(A.gset.cos.stat, 2),                      
                     ", cos.pvalA=",round(cos.pvalA, 4), 
                     ", dop.pvalB=", round(dop.pvalA, 4),
                     "; dop=", round(A.gset.dop.stat, 2),                     
                     ", cos.pvalB=", round(cos.pvalB,4), 
                     ", dop.pvalB=", round(dop.pvalB, 4))
    plot(x, y, xlab = x.lbl, ylab = y.lbl, pch = 19, col = "gray",
        sub = sub.lbl, main = gset.id)
    points(x[agr], y[agr], col = "red", pch = 19)
    lines(range(x), c(0, 0), lwd = 5, col = "black")
    lines(c(0, 0), range(y), lwd = 5, col = "black")
    pca.res<-prcomp(cbind(x,y),center=F)
    abline(a=0,b=pca.res$rotation[2,1]/pca.res$rotation[1,1],lwd=5,col="blue")    
       
    }
   
}

