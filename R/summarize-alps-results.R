rm(list=ls());

library(ape)
library(phytools)
library(geiger)
library(survival)
library(data.table)
library(stringr)

# set to results directory
wd <- '/Users/tomahern/Dropbox/ALPS/alps_121417/repository/1711-lash-tamoxifen-pathway-p1717801-master/'
setwd(wd)

header <- 'logs/results_simpath_cohort/'
header2 <- 'simpath-cohort-scram-coxph-ties-timer'

# read in raw output files

fstree <- fread(paste0(wd,header,header2,"-tree.txt"),sep=";",header=F)
fssampler <- fread(paste0(wd,header,header2,"-sampler.txt"),sep="\t",header=T)
fspsi <- fread(paste0(wd,header,header2,"-psi.txt"),sep="\t",header=T)
fsparameters <- fread(paste0(wd,header,header2,"-tree-parameters.txt"),sep="\t",header=T)
fsthetas <- fread(paste0(wd,header,header2,"-tree-thetas.txt"),sep="\t",header=T)

# Basic checks
# start by checking psi, the hyperparameter related to the prior forest 
print("Hyperparameter (psi)")
table(fspsi$psi)

# report the number of completed iterations
print("Number of iterations")
niter <- nrow(fssampler)
print(niter)

# check the acceptance rate for proposed pathway changes
acceptrate <- length(fssampler$accept[fssampler$accept==TRUE])/niter
print("Acceptance rate")
print(acceptrate)

# number of internal nodes
print("Number of internal nodes")
print(table(fssampler$nnodes)/niter)

# distance from prior forest
print("Distance from prior forest")
print(table(fssampler$distance))

# Summarize tree structures
# summarize the trees
# start by linking to sampler
fstree$V2 <- NULL
fstree$treeindex <- as.numeric(row.names(fstree))
fstree$tree <- substr(fstree$V1,regexpr('\\(',fstree$V1), nchar(fstree$V1))
fssamplertree <- merge(fssampler,fstree,by.x="treeindex",by.y="treeindex",all.x=T,sort=F)

# tree posteriors
treepost <- aggregate(fssamplertree$tree,by=list(fssamplertree$tree),FUN=function(x){NROW(x)})
names(treepost) <- c("tree","count")
treepost$post <- treepost$count/niter
treepost$perc <- treepost$post*100
treeprior <- mean(treepost$post)
treepost$BF <- round((treepost$post/(1-treepost$post))/(treeprior/(1-treeprior)),0)
print(treepost)

# Posterior odds for individual variants
# genetic factor marginal posteriors
text<-paste(treepost$tree,";",sep="")
temp <- read.tree(text=text)
getlabels <- function(tree) {
  return(tree$tip.label)
}
tips <- lapply(temp,getlabels)
nt <- unique(unlist(tips))
main <- matrix(0,ncol=length(nt))
names(main) <- nt
for (i in 1:nrow(treepost)) {
  main[tips[[i]]] <- main[tips[[i]]] + treepost$count[i]
}
mainpost <- data.frame(count=main,post=main/niter)

# calculate Bayes factors (BF)
mainprior <- mean(mainpost$post)
mainpost$bf <- round((mainpost$post/(1-mainpost$post))/(mainprior/(1-mainprior)),0)
print("Genetic factor marginal posteriors")
print(mainpost)

# pull out parameters for any tree with BF >= 1
treepost$perc <- round(treepost$post*100)
s <- subset(treepost,BF>=1)
print("Trees with BF >= 1")
print(s)

s2 <- subset(mainpost,bf>=1)

# plot cladograms of selected tree structures
par(mfrow=c(2,2))

for (i in 1:nrow(s)) {
  print(s$tree[i])
  ti <- min(fstree[fstree$tree==s$tree[i]]$treeindex)
  tiparam <- subset(fsparameters,treeindex==ti)
  print(tiparam)
  tithetas <- subset(fsthetas,treeindex==ti)
  print(tithetas)
  t <- read.tree(text=paste0(s$tree[i],";"))
  plot(t,root.edge=TRUE,main=paste0("Posterior: ",s$perc[i],"% BF:",s$BF[i]),type="cladogram",sub=paste0("beta=",round(tiparam$beta,2)))
  edgelabels(round(as.vector(t(subset(tithetas,select=c("theta.1","theta.2")))),2))
  #par(ask=TRUE) 
}

# Pathway concept ALPS summarization
# read in the spreadsheet with concept definitions
path <- read.csv('data/summarize_everything.csv',header=T)

# concept names start in column 3
path.names <- names(path)[3:ncol(path)]
for (i in path.names) {
  
  cur <- subset(path,!is.na(path[i]))
  
  # count the times any concept contributor appears in an ALPS tree
  concept.post <- 0
  for (j in 1:nrow(treepost)) {
    test <- any(tips[[j]]%in%cur$datavar)
    if(test) {
      concept.post <- concept.post + treepost$count[j]
    }	
  }
  concept.post <- concept.post/niter
  print(paste0("concept:", i , " posterior:", concept.post))
}
