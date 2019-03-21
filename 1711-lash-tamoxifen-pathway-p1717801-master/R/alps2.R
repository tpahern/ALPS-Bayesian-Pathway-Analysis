#-------------------------------------------------------------
# ALPS 2.0 Beta
# (c) 2017 James Baurley
#-------------------------------------------------------------

library(ape)
library(phytools)
library(geiger)
library(survival)
library(zoo)

#-------------------------------------------------------------
# Fit the model
#------------------------------------------------------------- 

fitalps <- function(iter,curtree,initipsi,normpotts,prioronly,lik,prefix) {
	# This function samples trees
	# Args:
 	#	iter: the number of MCMC iterations   
 	#	curtree: the initial tree
 	#	initipsi: where to start in the psi grid
 	#	normpotts: NA or TRUE, run without data to derive the normalizing constant for the potts prior
 	# 	prioronly: NA or TRUE, run without data to sample from the prior distribution
 	# 	prefix: adds to results filename
 	# 	lik: 'coxph-noties' for cox ph regression not accounting for ties, 'coxph-ties for breslow approximation (slower), 'logistic' for unconditional logistic regression
 	# TODO: Create seperate output file for beta's for logistic regression intercept and covariates
 	
	# Create mask for ties, not used unless there are ties
	mask <- array(dim=0)
	
	if (lik%in%c('coxph-ties','coxph-noties')) {
	
		# Assumes time.var is sorted!  Error if not
		try(if(is.unsorted(time.var)) stop("time var not sorted"))
		
		if(lik=='coxph-ties') {
			mask <- as.numeric(!duplicated(time.var))
			mask[mask==0] <- NA
		}
	}
	
	howlong <- iter;
	curtree$negloglik <- NA
	curtree$beta <- NA
	curtree$theta <- NA
	curtree$logpost <- NA
	
	curtree$psi <- psi[initipsi]
	curtree$ipsi <- initipsi
	
	usedata <- TRUE
	if (normpotts || prioronly) {
		usedata <- FALSE
	}
	
	# begin by fitting the starting tree
	if (usedata) {
		curtree <- optimtree(curtree,lik,mask)
	}
	curtree <- post(curtree,usedata)

	treeindex <- 1

	for (i in 1:howlong) {
		if (i%%10==0) print(paste0((i/howlong)*100,"% complete."))
		
		append <- TRUE
		accept <- FALSE
		loga <- NA
		ru <- NA
		if (i==1) append <- FALSE
	
		# propose a change
		proposedtree <- newtree(curtree)
		write.tree(proposedtree,file=paste0(wd,prefix,"-proposed-trees.txt"),append=append,tree.names=i)	
	
		proposedtree$psi <- curtree$psi
		proposedtree$ipsi <- curtree$ipsi
	
		# fit tree, compute posterior
		proposedtree$negloglik <- NA;
		proposedtree$beta <- NA;
		proposedtree$theta <- NA;
		proposedtree$logpost <- NA;
		
		if (usedata) {
			proposedtree <- optimtree(proposedtree,lik,mask)
		}
		proposedtree <- post(proposedtree,usedata)
		
		loga <- proposedtree$logpost - curtree$logpost
		loga <- loga + log(1/proposedtree$Q)
		ru <- log(runif(1))
		
		if (!is.na(loga)) {
			if (ru < loga) {
				curtree <- proposedtree;
				accept <- TRUE;
				treeindex <- treeindex+1;
			}
		}
		
		# update psi every iteration
		if (!normpotts) {
			curtree <- updatepsi(curtree,i,usedata,prefix)
		}
		
		# lots of logging
			
		df <- data.frame(i=i,treeindex=treeindex,accept=accept,nnodes=Nnode(curtree),distance=curtree$distance,loglik=-curtree$negloglik,post=curtree$logpost,loga=loga,Q=proposedtree$Q,ru=ru)
	
		write.table(df,file=paste0(wd,prefix,"-sampler.txt"),append=append,quote=F,sep="\t",row.names=F,col.names=!append)
			
		if (i==1 || accept==TRUE) {
			write.tree(curtree,file=paste0(wd,prefix,"-tree.txt"),append=append,tree.names=treeindex)	
			paramdf <- data.frame(treeindex=treeindex,nnodes=Nnode(curtree),beta=curtree$beta,negloglik=curtree$negloglik,logpost=curtree$logpost,distance=curtree$distance)
			write.table(paramdf,file=paste0(wd,prefix,"-tree-parameters.txt"),append=append,quote=F,sep="\t",row.names=F,col.names=!append)
			
			if (!anyNA(curtree$theta) && !is.null(curtree$theta) ) {
				thetadf <- data.frame(treeindex=treeindex,node=row.names(curtree$theta),theta=curtree$theta)
  				write.table(thetadf,file=paste0(wd,prefix,"-tree-thetas.txt"),append=append,quote=F,sep="\t",row.names=F,col.names=!append)
			}
		}			
	}
}

#-------------------------------------------------------------
# Computes the log posterior of a given tree
#-------------------------------------------------------------

post <- function(tree,usedata) {
	# This function computes the approximated posterior for a tree
	# Args:
 	#   tree: the current tree
 	#	usedata: TRUE if using data
  	# Returns:
  	#   tree: a tree with logpost and distance set
  	# TODO: re-consider how distance is calculated from the similarity score

	tree$logpost <- NA;
	tree$distance <- NA;
	#npar <- ncol(covariate.matrix)+2*Nnode(tree)  # covariates and 2 theta's for every node
	npar <- 2*Nnode(tree)  # covariates and 2 theta's for every node

	n <- ncol(dos)
	gamma <- 1.0 # for extended BIC, set to 0.5 or 1 for applications with many parameters

	# penalizes based on the number of model parameters and the number of trees possible for a given number of tips
	sj <- howmanytrees(length(tree$tip.label))
	#**** check if exists, if it doesn't set map to zero 
	map <- NA;
	
	# if not using data, set map to zero so can estimate prior
	if (!usedata) {
		map <- 0
	}
	
	if (!is.na(tree$negloglik)) {
		map <- -tree$negloglik
	}
	eBIC <- -2*map + npar*log(n) + 2*gamma*log(sj)

	potts <- 0;
	
	# creates a distance from the similarity score
	tree$distance <- 1/(priorsim(prior.forest,tree)+1)
	
	potts <- 2*tree$psi*tree$distance 
	
	logpost <- -0.5*(potts + eBIC);
	tree$logpost <- logpost;
	return(tree);
}

#-------------------------------------------------------------
# Generates a new tree from an existing one
#-------------------------------------------------------------

newtree <- function(tree) {
	# This generates a new tree from an existing one
	# Args:
 	#   tree: the current tree
  	# Returns:
  	#   tree: a modified tree with the Q-ratio set
  	# Todo: consider 1-input trees or prior on theta=c(0,1) or c(1,0) to prefer changes to topology rather than one factor having a small or null effect yet still in the tree.
  	 	
 	Q=1;
 	nodes <- Nnode(tree)
 	
 	# a 2-tip tree is a special case, no delete move
 	# a saturated tree is a special case, no add move, no nodes to replace, only option is delete
 	if (nodes==1) {
 		move <- sample(c('add','replace'),1)
 		Q <- Q * (1/2)/(1/3)
  	} else if (nodes==(nrow(dos)-1)) {
  		move <- 'delete'
  		Q <- Q * 1/(1/3)
  	} else {
 		move <- sample(c('add','delete','replace'),1)	
	}
		
	if (move=='add') {
 		tree = addtip(tree);
 		# moving to a saturated tree is a special case, you can only delete
 		pnodes <- Nnode(tree)
 		if (pnodes==(nrow(dos)-1)) {
 			Q <- Q * (1/3)/(1)
 		}
 	} else if (move=='delete') {
 		tree = deletetip(tree);
 		# moving back to a 2-tip tree is a special case, on the way back there is no delete
 		if (Nnode(tree)==1) {
 			Q <- Q * (1/3)/(1/2) 
 		}
 	} else  {
 		tree = replacetip(tree);
 	}
	#print(paste0("Q:",Q))
	
	tree$Q <- tree$Q * Q;
	return(tree);
}

#-------------------------------------------------------------
# Tree operations: addtip, deletetip, and replacetip
#-------------------------------------------------------------

addtip <- function(tree) {
	# This adds a tip to the tree
	# Args:
 	#   tree: the current tree
  	# Returns:
  	#   tree: a tree with a tip added and Q-ratio set
 
 	Q1 <- 1;  	#Pr(tree->proposedtree)
 	Q2 <- 1;	#Pr(proposedtree->tree)
 
 	# select where to add
  	len.tip <- length(tree$tip.label)
	which.tip <- sample(len.tip,1)
	Q1 <- Q1 * (1/len.tip)
	
	# pick variable to add
	# candidates are any variable except the ones currently in the tree
	candvars <- rownames(dos)
	candvars <- candvars[!(candvars%in%tree$tip.label)]
	newtip <- sample(candvars,1);
 	Q1 <- Q1 * (1/length(candvars))
 	
 	tree <- bind.tip(tree,tip.label=newtip,where=which.tip)

	Q2 <- Q2 * (1/(len.tip+1));

	#print(paste0("Add ", Q1, " ", Q2, " Qratio: ",(Q1/Q2)));
	tree$Q <- Q1/Q2    # set Q-ratio
	
	return(tree)
}

deletetip <- function(tree) {
  	# This deletes a tip from the tree
  	# Args:
  	#   tree: the current tree
  	# Returns:
  	#   tree: a tree with a tip deleted and the Q-ratio set
 	
 	Q1 <- 1;  	#Pr(tree->proposedtree)
 	Q2 <- 1;	#Pr(proposedtree->tree)
 
 	len.tip <- length(tree$tip.label)
	
	# for >2-tip tree, delete tip, otherwise return current tree with a warning
	 if (len.tip > 2) {
		which.tip <- sample(len.tip,1)
		tree <- drop.tip(tree,which.tip)		
	} else {
		warning("2-tip tree passed to deletetip")
	}
	Q1 <- Q1 * (1/len.tip)
	Q2 <- Q2 * (1/(len.tip-1)) # select where to add
	
	candvars <- rownames(dos)
	candvars <- candvars[!(candvars%in%tree$tip.label)]

	Q2 <- Q2 * (1/(length(candvars))); # pick a candidate variable to add
	
	#print(paste0("Add ", Q1, " ", Q2, " Qratio: ",(Q1/Q2)));
	
	tree$Q <- Q1/Q2    # set Q-ratio
	
	return(tree)
}

replacetip <- function(tree) {
  	# This replaces a tip with a new variable
  	# Args:
 	#   tree: the current tree
  	# Returns:
  	#   tree: a tree with the replaced tip and the Q-ratio set
 
 	# select tip to change
  	len.tip <- length(tree$tip.label)
	which.tip <- sample(len.tip,1)
	
	# pick variable to add
	# candidates are any variable except the ones currently in the tree
	candvars <- rownames(dos)
	candvars <- candvars[!(candvars%in%tree$tip.label)]
	newtip <- sample(candvars,1);

 	tree$tip.label[which.tip] <- newtip;
	
	tree$Q <- 1    # set Q-ratio
	
	return(tree)
}

#-------------------------------------------------------------
# Compare tree to prior forest
#-------------------------------------------------------------

priorsim <- function(forest,tree){
# This computes the similarity of the current tree to the prior forest
  	# Args:
 	#   forest: the user-specified trees
 	#	tree: the current tree
  	# Returns:
  	#   sim: the number of matching pairwise combinations
 
	sim <- 0
	tips <- tree$tip
	comb <- t(combn(tips,2))
	sim <- nrow(merge(comb,forest,by.x=c("V1","V2"),by.y=c("snp1","snp2")))
	return(sim)
}

                            
#-------------------------------------------------------------
# Tree-based model fitting
#-------------------------------------------------------------
 
alpstreeph <- function(pars, x, xn, y,lik,tiemask) {
	# This is the likelihood function for a ALPS tree, set up for the cox ph or logistic model
	# Args:
  	#   pars: The model parameters
 	#   x: the tree
 	#	xn: design matrix of observed (tips) and predicted values for the tree
  	#   y: the response variable
  	# 	lik: coxph-ties,coxph-noties or logistic
  	#	tiemask: for quick handling of ties, if applicable
  	# Returns:
  	#  the computed log likelihood
  	# TO DO: generalize for covariates
 
  	# creates a matrix of theta's from pars 
  	ntips <- length(x$tip);
	
  	theta <- matrix(pars[1:(2*Nnode(x))],ncol=2,byrow=T)
  	rownames(theta) <- (ntips+1):(Nnode(x)+ntips)  # labels the rows to match the tree

  	#beta1 <- pars[length(pars)] # net effect coefficient
  	beta <- pars[(2*Nnode(x)+1):length(pars)]
  	
 	# works upwards towards the tree root  
 	lastnode <- length(x$tip) +1 ;
 	start <- length(x$tip) + Nnode(x)

  	for (i in start:lastnode) {
  		# gets the children of the internal node i
  		edge <- x$edge[x$edge[,1]==i,]
  		i1 <- edge[1,2]
    	i2 <- edge[2,2]
  
      	tix <- as.character(i)
  		
  		# fill in the predicted values for each internal node
		xn[i,] <- theta[tix,1] * xn[i1, ] + theta[tix,2] * xn[i2, ] + (1 - theta[tix,1] - theta[tix,2]) * xn[i1, ] * xn[i2, ];  
  	}

	# the predicted values for the last node are the ones that are plugged into the regression model
	pred <- as.matrix(xn[lastnode,])
	
	
	# cox regression likelihood, assumes data is sorted
	if (lik%in%c('coxph-ties','coxph-noties')) {
		a <- pred %*% beta
		if (length(tiemask)==0) {
			b <- log(rev(cumsum(rev(exp(a)))))
		} else {
			b <- rev(cumsum(rev(exp(a))))
			b <- log(na.locf(b*tiemask))
		} #b2 <- ave(b, time.var, FUN = max)	
  		-sum((a - b)[y==1])
  	}
  	else if (lik%in%c('logistic')) {		 
		pred <- cbind(rep(1,nrow(pred)),pred) # add in intercept
		xb <- pred%*%beta
		#loglik <- sum(-y*log(1+exp(-xb)) - (1-y)*log(1+exp(xb)))
		p <- plogis(xb)
		p <- ifelse(p >= (1-.Machine$double.eps), (1-.Machine$double.eps), p)
		loglik <- sum(ifelse(y,log(p),log(1-p)))
		-loglik
    }
}  

optimtree <- function(tree,lik,mask) {
  	# This function optimizes the likelihood function (alpstree) for a ALPS tree.
  	#
  	# Args:
  	#   tree: the ALPS tree representing the selected risk factors
  	#	lik: the likelihood function to optimize
  	#
  	# Returns:
	#   tree: the fitted tree

	
	
 	# create design matrix xn for given tree of observed dosages (leaves) and predicted values (internal nodes)
	ntips <- length(tree$tip);
	xn <- matrix(nrow=(Nnode(tree)+ntips),ncol=ncol(dos))
	xn[1:ntips,] <- dos[tree$tip,] 
  
  	# create parameter vector for a given tree 
  	# the beta set of parameters are for the regression
  	# the theta set of parameters are for the system of equations, initialized them to 0.5
  	
 	theta <- array(0.5,dim=(2*Nnode(tree)));
    
    nbeta <- 1
    if (lik%in%c('logistic')) {
    		nbeta <- 2
    }
    
    beta <- array(0,dim=nbeta)
    if (lik%in%c('logistic')) {
    		beta[1] <- 0
    }

    pars <- c(theta,beta)
   	
   	# set lower constaints of zero for thetas, -inf for beta
  	lower <- array(-Inf,dim=length(pars))
  	lower[1:length(theta)] <- 0
  	
	# set upper constraints
	upper <- array(Inf,dim=length(pars))
	upper[1:length(theta)] <- Inf  # removing the constraint on theta
	
   	op <- try(optim(
   			pars,
      		fn = alpstreeph,
      		hessian = F,
      		method = "L-BFGS-B",
      		lower = lower,
      		upper = upper,
      		x = tree,
      		xn = xn,
      		y = status.var,
      		lik=lik,
      		tiemask=mask),silent=TRUE)   #; ,control=list(trace=6)
  	
	tree$beta = NA;
	tree$theta = NA;
	tree$negloglik = NA;
	
	if(class(op)=="try-error") {
		print("WARNING: Finite value error in optimizing tree...moving on")
	} else if (op$convergence == 0) {
  		tree$beta <- op$par[length(pars)]	# pathway effect is always the last element of par
  		tree$theta <- matrix(op$par[1:length(theta)],ncol=2,byrow=T)
		rownames(tree$theta) <- (ntips+1):(Nnode(tree)+ntips)
		tree$negloglik <- op$value
  	}

  	 	return(tree)
}



#---------------------------------------------------------------------------------
# Updates the psi parameter for how quickly the prior fades as distance increases
#---------------------------------------------------------------------------------

updatepsi <- function(tree,i,usedata,prefix) {
 	# This function does a random walk on the psi parameter in the potts prior
  	#
  	# Args:
  	#   tree: the current tree
  	#	i: the current MCMC iteration
  	#	usedata: indicator if data is being used in this run	
  	# 	prefix: file prefix for writing out psi
  	# Returns:
	#   tree: the tree with psi parameter updated
	#TODO: move write out of psi to sampler outputs

	# if index is on either end of the grid move
	ipsi <- tree$ipsi
	
	if (ipsi==1) {
		newipsi <- 2
		qratio <- 2
	} else if (ipsi==npsi) {
		newipsi <- npsi-1;
		qratio <- 2
	} else if (runif(1) < 0.5) {   # otherwise move randomly
		newipsi <- ipsi-1
		# if not moving to an end
		if (ipsi>2 && ipsi<(npsi-1)) {
			qratio <- 1
		} else {
			qratio <- 0.5
		}
	} else {
		newipsi <- ipsi+1
		if (ipsi>2 && ipsi<(npsi-1)) {
			qratio <- 1
		} else {
			qratio <- 0.5
		}
	} 
	
	# compute old potts prior dividing by normalizing constant
	oldpotts <- 2*tree$psi*tree$distance - 2*log(cpsi[ipsi])
	
	# compute new psi and potts prior 
	newpsi <- psi[newipsi]
	newpotts <- 2*newpsi*tree$distance - 2*log(cpsi[newipsi])
	
	# hasting step
	likratio <- newpotts - oldpotts
	priorratio <- log(tree$psi/newpsi)
	hratio <- likratio + priorratio + log(1/qratio) 
	accept=FALSE
	# if accepted
	if (hratio < log(runif(1))) {
		accept=TRUE
		tree$psi <- newpsi
		tree$ipsi <- newipsi
		tree <- post(tree,usedata) # maybe look for more efficient way (recomputes distance)
	}
	append <- TRUE
	if (i==1) append <- FALSE
	
	psidf <- data.frame(i=i,distance=tree$distance,psi=tree$psi,newpsi=newpsi,ipsi=ipsi,newipsi=newipsi,likratio=likratio,priorratio=priorratio,hratio=hratio,accept=accept)
	
	write.table(psidf,file=paste0(wd,prefix,"-psi.txt"),append=append,quote=F,sep="\t",row.names=F,col.names=!append)
	
	return(tree)
}

