#-------------------------------------------------------------
# Simulated cohort ALPS analysis
#-------------------------------------------------------------

# Load needed libraries
library(ape)
library(phytools)
library(geiger)
library(survival);

rm(list=ls());


# sets the seed so we can reproduce the run
set.seed(1);

# load in alps
source(file=paste0(wd,"alps2-test-optim.R"))

# Read in Dataset
ds <- read.csv(paste0(wd,'cacosim.csv'),header=T)

# Read in SNP-level prior specification
all <- read.table(file=paste0(wd,"prior-forest-snps.txt"),sep="\t",header=T)
prior.forest <- subset(all,select=c("snp1","snp2"));

# Config dataset for ALPS
#time.var <- ds$time_tie
#status.var <- ds$event
case <- ds$case

# Note that it is called dos for genotype dosages. 
dos <- ds;
dos$case <- NULL
dos$strat <- NULL
#dos$event <- NULL

dos <- t(as.matrix(dos))

#-------------------------------------------------------------
# Initialize ALPS and run
#-------------------------------------------------------------
load(file=paste0(wd,'cpsi.RData'))

# start at a random spot in the prior forest

spot <- prior.forest[sample(nrow(prior.forest),1),]
spot.tree <- paste0("(",spot$snp1,",",spot$snp2,");")
curtree <- read.tree(text=spot.tree)

# Tom, change iterations.  I would see how long 1000 takes, then scale up. Goal is to get 100,000.
#system.time(fitalps(iter=100000,initipsi=19,curtree,normpotts=FALSE,prioronly=FALSE))


### try to write optim function
#x <- dos[12,]
x <- dos[12,]

rfit <- glm(case ~ x,family=binomial)
logLik(rfit)

logisticloglik <- function(beta,xn,y) {
	pred <- as.matrix(xn) 
	pred <- cbind(rep(1,nrow(pred)),pred)
	xb <- pred%*%beta
	loglik <- sum(-y*log(1+exp(-xb)) - (1-y)*log(1+exp(xb)))
	return(-loglik)	
}


system.time(myfit <- optim(par=c(0,0),fn=logisticloglik,xn=x,y=case,method="L-BFGS-B",lower=c(-Inf,-Inf),upper=c(Inf,Inf)))

#logisticloglik(2,x,case)






# old stuff

system.time(fit <- coxph(Surv(time.var,status.var) ~ x,ties="breslow"))
logLik(fit)

survdf <- data.frame(time.var=time.var,status.var=status.var,x=x)
# order data by increasing time
survdf <- survdf[order(survdf$time.var),]

likfun<- function(beta, data) {
  X <- as.matrix(data[, -(1:2)])
  a <- X %*% beta
  b <- log(rev(cumsum(rev(exp(a)))))
b2 <- with(data, ave(b, time.var, FUN = max))

  -sum((a - b2)[data$status.var==1])
}


 #t <- as.matrix(unique(survdf$time.var[survdf$status.var==1]))
 t <- unique(survdf$time.var[survdf$status.var==1])
 
 event <-  subset(survdf,status.var==1)

likbres <-function(beta) {
	# breslow approx
	ll <- 0
	for (i in 1:length(t)) {
 		Dm <- subset(event,time.var==t[i])
 		Rm <- subset(survdf,time.var>=t[i])

 		Dmb <- sum(Dm$x*beta)
 		eRmb <- sum(exp(Rm$x*beta))
 		#eDmb <- sum(exp(Dm$x*beta))
 	
 		#ks <- 0:(nrow(Dm)-1)/nrow(Dm)
 
 		partb <- nrow(Dm)*log(eRmb)
 	
 		ll <- ll + Dmb - partb
}
 return(-ll)
}
	
logLik(fit)
system.time(temp<-nlm(f=likbres, p=0))
system.time(temp2<-nlm(f=likfun,p=0,data=survdf))





likfun2 <- function(beta) {
# efron approx
 ll <- 0;
 for (i in 1:length(t)) {
 	Dm <- subset(event,time.var==t[i])
 	Rm <- subset(survdf,time.var>=t[i])

 	Dmb <- sum(Dm$x*beta)
 	eRmb <- sum(exp(Rm$x*beta))
 	eDmb <- sum(exp(Dm$x*beta))
 	
 	ks <- 0:(nrow(Dm)-1)/nrow(Dm)
 
 	partb <- sum(log(eRmb-ks*eDmb))
 	
 	ll <- ll + Dmb - partb
}
 return(-ll)
}

likfun3 <- function(beta){
	blaa <- unlist(lapply(talt,FUN=appl,beta=beta))
	return(-sum(blaa))
}


appl <- function(tr,beta) {
	Dm <- event[event$time.var==tr,]
 	Rm <- survdf[survdf$time.var>=tr,]
    Dnbpart <- Dm$x*beta
 	Dmb <- sum(Dnbpart)
 	eRmb <- sum(exp(Rm$x*beta))
 	eDmb <- sum(exp(Dnbpart))
 	
 	ks <- 0:(nrow(Dm)-1)/nrow(Dm)
 
 	partb <- sum(log(eRmb-ks*eDmb))
 	return(Dmb-partb)
}




system.time(temp<-nlm(f=likfun3, p=0))
system.time(temp2<-nlm(f=likfun2, p=0))


temp <- nlm(f=likfun,p=1,data=survdf)

library(optimx)
system.time(temp2 <- optim(par=0,fn=likfun3,method="L-BFGS-B",lower=-Inf,upper=Inf))

function(beta) {
# efron approx
 ll <- 0;
 for (i in 1:length(t)) {
 	Dm <- subset(event,time.var==t[i])
 	Rm <- subset(survdf,time.var>=t[i])
 	parta <- sum(Dm$x*beta)
 	partb <- 0
 	for (k in 0:(nrow(Dm)-1)) {
 		a <- sum(exp(Rm$x*beta))
 		b <- (k/nrow(Dm)) * sum(exp(Dm$x*beta))
 	 	partb <- partb+log(a-b)
 	}
 	ll <- ll + parta - partb
}
 return(-ll)
}

# Conclusion: ALPS implementation in R, optimization of efron likelihood is not efficient. Could explore callout to C++ or C later

# Will explore with beslow approach





