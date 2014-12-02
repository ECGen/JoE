###JoE manuscript analysis
###MKLau 2may2014
library(sna)
library('RColorBrewer')
source('../src/func_town.R')
                                        #data pre-processing
r <- joePpd(read.csv('../data/joe_pwr_raw.csv'))[c(5,2,3,4,1,6,7,8,9),c(5,2,3,4,1,6,7,8,9)]
rs <- joePpd(read.csv('../data/joe_pwr_raw_sig.csv'))[c(5,2,3,4,1,6,7,8,9),c(5,2,3,4,1,6,7,8,9)]
m <- joePpd(read.csv('../data/joe_pwr_max.csv'))[c(5,2,3,4,1,6,7,8,9),c(5,2,3,4,1,6,7,8,9)]
ms <- joePpd(read.csv('../data/joe_pwr_max_sig.csv'))[c(5,2,3,4,1,6,7,8,9),c(5,2,3,4,1,6,7,8,9)]
r <- r[rownames(r)!='TreeWebs.09',colnames(r)!='TreeWebs.09']
rs <- rs[rownames(rs)!='TreeWebs.09',colnames(rs)!='TreeWebs.09']
m <- m[rownames(m)!='TreeWebs.09',colnames(m)!='TreeWebs.09']
ms <- ms[rownames(ms)!='TreeWebs.09',colnames(ms)!='TreeWebs.09']
                                        #plot
v.col <- c(grey(0.25),rep("#662506",3),rep("#254C00",4))

new.coord <- FALSE
if (new.coord){coord <- locator(nrow(r));dput(coord,file='../data/coords.rdata')}else{
coord <- dget(file='../data/coords.rdata')}
coord$y[1] <- 1.27281660
coord$x[7] <- -1.45703700
labels <- c('Epiphytic Lichens 2010','Ectomycorrhizal Fungi 2006','Soil Bacteria 2004','Soil Fungi 2004','Twig Endophytes 2006','Leaf/Stem Modifying Arthropods 2010','Pathogens 2009','Pathogens 2010')
###
pdf(file='~/projects/JoE/results/netFig2a.pdf',width=10,height=6.5)
par(family='Times',mai=rep(0,4))
joePlot(abs(r),rs,v.col=v.col,v.cex=1)
text(x=c(-1.457897,-1.034757),y=c(0.0765118,1.084242),labels='Fungal Leaf',cex=1.25)
text(x=coord$x,y=coord$y,labels=labels,cex=1.25)
legend('topleft',legend='A',bty='n',cex=2.5)
dev.off()
pdf(file='~/projects/JoE/results/netFig2b.pdf',width=10,height=6.5)
par(family='Times',mai=rep(0,4))
joePlot(abs(m),ms,v.col=v.col,v.cex=1)
text(x=c(-1.457897,-1.034757),y=c(0.0765118,1.084242),labels='Fungal Leaf',cex=1.25)
text(x=coord$x,y=coord$y,labels=labels,cex=1.25)
legend('topleft',legend='B',bty='n',cex=2.5)
dev.off()
pdf(file='~/projects/JoE/results/netFig2L.pdf')
par(family='Times',mai=rep(0,4))
plot(1:10,pch='',xaxt='n',yaxt='n',bty='n',xlab='',ylab='')
legend('center',legend=c(expression(italic('P')<=0.05),expression(paste('0.05 < ',italic('P'),' < ','0.1')),expression(italic('P')>=0.1),'','Phyllosphere','Trunk','Soil'),pch=c(15,15,15,15,19,19,19),col=c('black',grey(0.5),grey(0.85),0,'#254C00',grey(0.25),"#662506"),cex=3,box.lwd=1,box.col='grey')
dev.off()

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
