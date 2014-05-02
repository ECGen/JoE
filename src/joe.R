###JoE manuscript analysis
###MKLau 2may2014

library(sna)
library(enaR)
                                        #
x <- read.csv('../data/mantelRhoTable.csv')
rownames(x) <- x[,1]
x <- x[,-1]
x <- rbind(rep(0,ncol(x)),x)
x <- cbind(x,rep(0,nrow(x)))
rownames(x)[1] <- colnames(x)[1]
colnames(x)[ncol(x)] <- rownames(x)[nrow(x)]
x[is.na(x)] <- 0
x <- as.matrix(x)
x[upper.tri(x)] <- t(x)[upper.tri(x)]
x.val <- x
                                        #
x <- read.csv('../data/mantelRhoSig.csv')
rownames(x) <- x[,1]
x <- x[,-1]
x[x<99] <- 0
x <- rbind(rep(0,ncol(x)),x)
x <- cbind(x,rep(0,nrow(x)))
rownames(x)[1] <- colnames(x)[1]
colnames(x)[ncol(x)] <- rownames(x)[nrow(x)]
x[is.na(x)] <- 0
x <- as.matrix(x)
x[upper.tri(x)] <- t(x)[upper.tri(x)]
x[x==99] <- -1
x[x==999] <- 1
x.sig <- x
                                        #
x <- abs(x.val)
x[abs(x.sig)!=1] <- 0
v.col <- c("green","red","red","red","green","green","green","green","green")
v.cex <- degree(x)
e.col <- x.sig
e.col[x.sig==-1] <- 'lightgrey'
                                        #

gplot(x,gmode='graph',displaylabels=TRUE,edge.lwd=(x*10)^0.85,vertex.col=v.col,vertex.cex=v.cex,vertex.border='white',edge.col=e.col)
                                        #without foliar arthropods
gplot(x[-9,-9],gmode='graph',displaylabels=TRUE,edge.lwd=(abs(x[-9.-9])*10)^0.85,vertex.col=v.col[-9],vertex.cex=v.cex[-9],vertex.border='white',edge.col=e.col)
