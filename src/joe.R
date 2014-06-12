###JoE manuscript analysis
###MKLau 2may2014
library(sna)
                                        #data pre-processing
r <- joePpd(read.csv('../data/joe_pwr_raw.csv'))[c(5,2,3,4,1,6,7,8,9),c(5,2,3,4,1,6,7,8,9)]
rs <- joePpd(read.csv('../data/joe_pwr_raw_sig.csv'))[c(5,2,3,4,1,6,7,8,9),c(5,2,3,4,1,6,7,8,9)]
m <- joePpd(read.csv('../data/joe_pwr_max.csv'))[c(5,2,3,4,1,6,7,8,9),c(5,2,3,4,1,6,7,8,9)]
ms <- joePpd(read.csv('../data/joe_pwr_max_sig.csv'))[c(5,2,3,4,1,6,7,8,9),c(5,2,3,4,1,6,7,8,9)]
                                        #plot
v.col <- c("purple","red","red","red","green","green","green","green","purple")

my.coord <- joePlot(r,rs,v.col=v.col,v.cex=2)
joePlot(m,ms,coord=my.coord,v.col,v.cex=2)
legend('topright',legend=c('P<0.05','0.05<P<0.10','P>0.05','','negative'),lty=c(1,1,1,1,2),lwd=c(10,3,1,0,1))

###Analytics for network
r. <- r;r.[rs!=2] <- 0;m. <- m;m.[ms!=2] <- 0
r.. <- r.[upper.tri(r.)][r.[upper.tri(r.)]!=0|m.[upper.tri(m.)]!=0]
m.. <- m.[upper.tri(m.)][r.[upper.tri(r.)]!=0|m.[upper.tri(m.)]!=0]
hist((r-m)[upper.tri(r)]);t.test((r-m)[upper.tri(r)]);wilcox.test((r-m)[upper.tri(r)])
hist((r.-m.)[upper.tri(r.)]);t.test((r.-m.)[upper.tri(r.)]);wilcox.test((r.-m.)[upper.tri(r.)])
hist((r..-m..));t.test((r..-m..));wilcox.test((r..-m..))

####
n <- x
n[x.sig<=0] <- 0
n <- n[-9,-9]
ev.cen <- round(evcent(n),3)
de.cen <- degree(n,'freeman')
names(ev.cen) <- names(de.cen) <- colnames(n)
par(mfrow=c(1,2))
barplot(ev.cen,las=2)
barplot(de.cen,las=2)

###Pathway proliferation
###At what path length is everything connected
##The largest strongly connected component
my.lcc <- function(x){
  x[eigen(sign(x))$vectors[,1]!=0,eigen(sign(x))$vectors[,1]!=0]
}

my.scc <- function(x){
  x[eigen(sign(x))$vectors[,1]==0,eigen(sign(x))$vectors[,1]==0]
}

n.lcc <- my.lcc(n)

pp <- sign(n.lcc)
i <- 2
while(any(pp==0)){
  pp <- pp + pp%^%i
  i <- i + 1
}

###
npp <- n.lcc%^%2
for (i in 3:4){
  npp <- npp * npp%^%i
}

###Multiplicative paths
mpp <- n*0

for (i in 1:nrow(mpp)){
  for (k in 1:nrow(mpp)){
    for (j in 1:nrow(mpp)){
      mpp[i,k] <- mpp[i,k]+ (n[i,j] * n[j,k])
    }
  }
}

gplot(mpp,gmode='graph',displaylabels=TRUE,vertex.col=v.col,vertex.border='white')

###If we add in the tree
a <- sign(n)
a <- cbind(a,c(rep(1,8)))
a <- rbind(a,c(rep(1,8),0))
rownames(a)[9] <- colnames(a)[9] <- 'Tree'


###Look at indirect path correlations 
library(expm)
out <- list()
for (i in 1:100){
  out[[i]] <- n%^%i
}
test <- lapply(out,function(x) x[1,2])
