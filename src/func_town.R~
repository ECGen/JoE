### Functions for Lamit
### Journal of Ecology invited paper analyses
### 17 Mar 2014
### MK Lau

###Network figs
joeLabels <- function(coord,labels){
  labels
  
}

joePpd <- function(x){
                                        #data pre-processing
  rownames(x) <- x[,1]
  x <- x[,-1]
  x[upper.tri(x)] <- t(x)[upper.tri(x)]
  diag(x) <- 0
  as.matrix(x)
}

joePlot <- function(x,y,coord,v.col,v.cex){
                                        #plots
  if (missing(v.cex)){v.cex <- 1}
  if (missing(v.col)){v.col <- rep('grey',nrow(x))}
  x. <- x.. <- x... <- x
  x.[y!=2] <- 0    
  x..[y!=1] <- 0
  x...[y!=0] <- 0
  e.lwd <- sqrt(abs(x)*100)
  e.col <- y;e.col[y==2] <- 'white';e.col[y==1] <- 'grey';e.col[y==0&x!=0] <- grey(0.95);e.col[e.col=="0"] <- "white"
  set.seed(0)
  coord <- gplot((x...),gmode='graph',displaylabels=FALSE,edge.lwd=e.lwd,
        vertex.col=0,vertex.cex=v.cex,vertex.border='white',
        label.bg='white',boxed.labels=TRUE,
        label.border=0,mode='circle',edge.col=e.col,
        edge.lty=1,edge.lty.n=2)
  set.seed(0)
  gplot((x..),gmode='graph',displaylabels=FALSE,edge.lwd=e.lwd,
        vertex.col=v.col,vertex.cex=v.cex,vertex.border='white',
        label.bg='white',boxed.labels=TRUE,
        label.border=0,mode='circle',
        edge.col='grey',new=FALSE,coord=coord,
        edge.lty=1,edge.lty.n=2)
  set.seed(0)
  gplot((x.),gmode='graph',displaylabels=FALSE,edge.lwd=e.lwd,
        vertex.col=v.col,vertex.cex=v.cex*2,vertex.border='white',
        label.bg='white',boxed.labels=TRUE,
        label.border=0,mode='circle',
        edge.col=1,new=FALSE,coord=coord,
        edge.lty=1,edge.lty.n=2)
  coord
}

### matrix = community matrices
### factor = genotype vector
### names = tree names
### matchComs: returns a list containing the 
### matrices and factors for the two inputs
### "matched" so that both only contain observations
### replicated a minimum number of times (DEFAULT = 2).
### 
### nMismatch: returns a matrix containing the number
### of names that do match across two sets of factors.
### xy = the number in x that are not in y (for a given level)
### yx = the number in y that are not in x (for a given level)
                                        #match communities based on a minimum number of reps
matchComs <- function(x='matrix',y='matrix',fx='factor x',fy='factor y',min=2){
  cx <- sapply(fx,function(x,y) table(y)[names(table(y))==x],fx)
  cy <- sapply(fy,function(x,y) table(y)[names(table(y))==x],fy)
  x <- x[cx>min,];y <- y[cy>min,]
  fx <- fx[cx>min];fy <- fy[cy>min]
  x <- x[fx%in%fy,];y <- y[fy%in%fx,]
  fx <- fx[fx%in%fy];fy <- fy[fy%in%fx]
  return(list(x=x,y=y,fx=fx,fy=fy))
}
                                        #
nMismatch <- function(x='names',y='names',fx='factor x',fy='factor y'){
  if (length(unique(fx))!=length(unique(fy))){warning('Factor levels not matched!')}
  fx <- as.character(fx);fy <- as.character(fy)
  xy <- numeric(length=length(unique(c(fx,fy))))
  yx <- numeric(length=length(unique(c(fx,fy))))
  for (i in 1:length(unique(fx))){
    xy[i] <- sum(sign((x[fx==unique(fx)[i]]%in%y[fy==unique(fy)[i]]==FALSE)))
    yx[i] <- sum(sign((y[fy==unique(fy)[i]]%in%x[fx==unique(fx)[i]]==FALSE)))
  }
  out <- rbind(xy,yx)
  colnames(out) <- as.character(unique(c(fx,fy)))
  return(out)
}

## x = tree names for dataset 1
## y = tree names for dataset 2
## fx = genotype list for dataset 1
## fy = genotype list for dataset 2

x <- factor(c('N1.1','N1.12','N2.13','N7.2'))
y <- factor(c('N1.1','N3.5','N2.13','N7.2'))
fx <- factor(c('T15','T15','RM2','RM2'))
fy <- factor(c('T15','T15','RM2','RM2'))
nMismatch(x,y,fx,fy)

###matchMats = match two matrices using 
## x = matrix
## y = matrix
## match.col = column number containing the matching information (DEFAULT = 1)

matchMats <- function(x,y,match.col=1){
  xiny <- sapply(x[,match.col],function(x,y) any(x%in%y),y=y[,match.col])
  yinx <- sapply(y[,match.col],function(x,y) any(x%in%y),y=x[,match.col])
  return(list(x[xiny,],y[yinx,]))
}
