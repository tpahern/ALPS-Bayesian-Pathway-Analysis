#-------------------------------------------------------------
# ALPS 2.0 Development
# (c) 2017 James Baurley
#-------------------------------------------------------------

library(ape)
library(phytools)
library(geiger)
library(survival);

#-------------------------------------------------------------
# Fit the model
#------------------------------------------------------------- 

fitalps <- function(iter,curtree,initipsi,normpotts,prioronly) {
	# This function samples trees
	# Args:
 	#	iter: the number of MCMC iterations   
 	#	curtree: the initial tree
 	#	normpotts NA or TRUE
 	# 	prioronly: NA or TRUE
 
	howlong <- iter;
	curtree$negloglik <- NA;
	curtree$beta <- NA;
	
	# start in middle for psi
	#ipsi <- round(length(cpsi)/2)
	# also tried low and high end and seems to wander appropriately
	
	curtree$psi <- psi[initipsi]
	curtree$ipsi <- initipsi
	
	# begin by fitting the starting tree
	if (!normpotts & !prioronly) {
		curtree <- optimtree(curtree);
	}
	curtree <- post(curtree);

	treeindex <- 1;

	for (i in 1:howlong) {
		if (i%%10==0) print(paste0((i/howlong)*100,"% complete."))
	
		append <- TRUE
		accept <- FALSE
		loga <- NA
		ru <- NA
		if (i==1) append <- FALSE
	
		# propose a change
		proposedtree <- newtree(curtree)
		proposedtree$psi <- curtree$psi
		proposedtree$ipsi <- curtree$ipsi
		
		# fit tree, compute posterior
		proposedtree$negloglik <- NA;
		proposedtree$beta <- NA;
		if (!normpotts & !prioronly) {
			proposedtree <- optimtree(proposedtree)
		}
		proposedtree <- post(proposedtree)

		loga <- proposedtree$logpost - curtree$logpost
		loga <- loga + log(1/proposedtree$Q)
		ru <- log(runif(1))
		
		if (ru < loga) {
			curtree <- proposedtree;
			accept <- TRUE;
			treeindex <- treeindex+1;
		}
		
		# update psi every iteration
		if (!normpotts) {
			curtree <- updatepsi(curtree,i)
		}
		
		# lots of logging, add psi here
	
		df <- data.frame(i=i,treeindex=treeindex,accept=accept,nnodes=Nnode(curtree),distance=curtree$distance,loglik=-curtree$negloglik,post=curtree$logpost,loga=loga,Q=proposedtree$Q,ru=ru)
	
		write.table(df,file=paste0(wd,"/results/testing-sampler.txt"),append=append,quote=F,sep="\t",row.names=F,col.names=!append)

		if (i==1 || accept==TRUE) {
			write.tree(curtree,file=paste0(wd,"/results/testing-sampler-tree.txt"),append=append,tree.names=treeindex)	
			paramdf <- data.frame(treeindex=treeindex,nnodes=Nnode(curtree),beta=curtree$beta,negloglik=curtree$negloglik,logpost=curtree$logpost,distance=curtree$distance)
			write.table(paramdf,file=paste0(wd,"/results/testing-sampler-tree-parameters.txt"),append=append,quote=F,sep="\t",row.names=F,col.names=!append)
			
			if (!is.na(curtree$theta)) {
				thetadf <- data.frame(treeindex=treeindex,node=row.names(curtree$theta),theta=curtree$theta)
  				write.table(thetadf,file=paste0(wd,"/results/testing-sampler-tree-thetas.txt"),append=append,quote=F,sep="\t",row.names=F,col.names=!append)
			}
		}			
	}
}

#-------------------------------------------------------------
# Computes the log posterior of a given tree
# 	Todo: revisit eBIC penalty
#	Todo: revisit Potts prior based on distance 
#		  (particularly fixed psi)
#-------------------------------------------------------------

post <- function(tree) {
	# This function computes the approximated posterior for a tree
	# Args:
 	#   tree: the current tree
  	# Returns:
  	#   tree: a tree with logpost and distance set

	tree$logpost <- NA;
	tree$distance <- NA;
	#npar <- ncol(covariate.matrix)+2*Nnode(tree)  # covariates and 2 theta's for every node
	npar <- 2*Nnode(tree)  # covariates and 2 theta's for every node

	n <- ncol(dos)
	gamma <- 0 # for extended BIC, set to 0.5 or 1 for applications with many parameters

	# penalizes based on the number of model parameters and the number of trees possible for a given number of tips
	sj <- howmanytrees(length(tree$tip.label))
	#**** check if exists, if it doesn't set map to zero 
	map <- 0;
	if (!is.na(tree$negloglik)) {
		map <- -tree$negloglik
	}
	eBIC <- -2*map + npar*log(n) + 2*gamma*log(sj)

	potts <- 0;
	tree$distance <- 1/(priorsim(prior.forest,tree)+1)
	
	potts <- 2*tree$psi*tree$distance 

	#if(exists('prior.forest')) {
	#tree$distance <- priordistance(prior.forest,tree)
	#	potts <- 2*psi*distance # check Potts paper 
	#}
	
	logpost <- -0.5*(potts + eBIC);
	tree$logpost <- logpost;
	return(tree);
}

#-------------------------------------------------------------
# Generates a new tree from an existing one
# 	Todo: add in PEAK probabilities to sample function
#-------------------------------------------------------------

newtree <- function(tree) {
	# This generates a new tree from an existing one
	# Args:
 	#   tree: the current tree
  	# Returns:
  	#   tree: a modified tree with the Q-ratio set
 	
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
 	Q2 <- 1;		#Pr(proposedtree->tree)
 
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
 	Q2 <- 1;		#Pr(proposedtree->tree)
 
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
	sim <- 0
	tips <- tree$tip
	comb <- t(combn(tips,2))
	sim <- nrow(merge(comb,forest,by.x=c("V1","V2"),by.y=c("snp1","snp2")))
	return(sim)
}


priorsim2 <- function(forest,tree) {
	sim <- 0
	for (i in 1: length(forest)) {
		matchedge <- which.edge(tree,forest[[i]]$tip)
		len <- length(matchedge)
			if (len > 1) {
				sim <- sim + 1
				#print(paste0(i," ", ptrees[[i]]$tip.label))
			}
	}
	return(sim)
}

priordistance <- function(forest, tree) {
 	# Computes the distance of a tree to a forest of priors specified by the user
  	# Args:
 	#	  forest: the collection of trees representing prior knowledge
 	#   tree: the current tree
 	# Returns:
  	#   dist: the computed distance
  	
  	# Todo: re-visit computation of distance
 
	gene.tree <- tree;
	tip.labels <- gene.tree$tip.label
	gene.tree$tip.label <- var[tip.labels,]$annot.num;

	#plot(gene.tree)
	#prior <- read.tree(text="((2,4),5);")

	# count matching edges between current tree and the prior forest
	sim <- 0;
	for (i in 1: length(forest)) {
		matchedge <- which.edge(gene.tree,forest[[i]]$tip)
		sim <- sim + length(matchedge)
	}
	# normalize similarity to the size of the current tree
	max <- length(forest)*nrow(gene.tree$edge)
	nsim <- sim/max;
	# compute distance
	dist <- 1-nsim;

	#print(paste0("sim:",sim," max:",max," nsim:",nsim, "dist:",dist))

	### plots
	#a <- which.edge(cur.gene.tree,prior.forest[[11]]$tip)
	#clcolr <- rep("darkgrey",dim(cur.gene.tree$edge)[1])
	#clcolr[a] <- "black"
	#plot(cur.gene.tree,lwd=3,edge.color=clcolr)

	return(dist)	
}


#-------------------------------------------------------------
# Tree-based model fitting
#-------------------------------------------------------------
 
alpstree <- function(pars, x, xn, y) {
	# This is the likelihood function for a ALPS tree, set up for an identity link function (linear regression).
  	# Args:
  	#   pars: The model parameters
 	#   x: the tree
 	#	xn: design matrix of observed (tips) and predicted values for the tree
  	#   y: the response variable
  	# Returns:
  	#   logl: the computed log likelihood
 
  	# creates a matrix of theta's from pars 
  	ntips <- length(x$tip);
	
  	theta <- matrix(pars,ncol=2,byrow=T)
  	rownames(theta) <- (ntips+1):(Nnode(x)+ntips)  # labels the rows to match the tree
  
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
	#pred <- cbind(1,xn[lastnode,])
	pred <- xn[lastnode,]
	
	# cox regression for matched case control
	# TODO: explore ways to speed this up 
	#		   	(1) Alternative [R] packages for study design)
	# 			(2) Handwrite likelihood function
	#			(3) Cache visited trees
	
  	reg <- coxph(Surv(time.var,status.var) ~ pred)

	negloglik <- -logLik(reg)[[1]]
 	beta1 <<- reg$coef[["pred"]]
  
  	# returns the -loglikelihood
  	return(negloglik)
}                            

optimtree <- function(tree) {
  	# This function optimizes the likelihood function (alpstree) for a ALPS tree.
  	#
  	# Args:
  	#   tree: the ALPS tree representing the selected risk factors
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
  	
  #  beta1 <- NA #pathway effect
  	theta <- array(0.5,dim=(2*Nnode(tree)));
  	
   	# set lower constaints of zero for thetas
  	lower <- array(-Inf,dim=length(theta))
  	lower[1:length(theta)] <- 0;
  	
	# set upper constraints  
	upper <- array(Inf,dim=length(theta)) 
	
   	op <- optim(
   		theta,
      	fn = alpstree,
      	hessian = F,
      	method = "L-BFGS-B",
      	lower = lower,
      	upper = upper,
      	x = tree,
      	xn = xn,
      	y = status.var);
  

	tree$beta = NA;
	tree$theta = NA;
	tree$negloglik = NA;

  	if (op$convergence == 0) {
  		tree$beta <- get("beta1");
  		tree$theta <- matrix(op$par,ncol=2,byrow=T)
		rownames(tree$theta) <- (ntips+1):(Nnode(tree)+ntips)
		tree$negloglik <- op$value
  	}
 	return(tree)
}

#---------------------------------------------------------------------------------
# Updates the psi parameter for how quickly the prior fades as distance increases
#---------------------------------------------------------------------------------

updatepsi <- function(tree,i) {
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
		tree <- post(tree) # maybe look for more efficient way (recomputes distance)
	}
	append <- TRUE
	if (i==1) append <- FALSE
	
	#print(paste0("i: ",i," distance: ",tree$distance, " psi: ", tree$psi, " newpsi: ",newpsi, " ipsi: ",ipsi," newipsi:", newipsi, " likratio:",likratio, " priorratio:", priorratio," hratio: ",hratio, " accept: ", accept))
	
	psidf <- data.frame(i=i,distance=tree$distance,psi=tree$psi,newpsi=newpsi,ipsi=ipsi,newipsi=newipsi,likratio=likratio,priorratio=priorratio,hratio=hratio,accept=accept)
	
	
	write.table(psidf,file=paste0(wd,"/results/testing-sampler-psi.txt"),append=append,quote=F,sep="\t",row.names=F,col.names=!append)
	
	return(tree)
	
}

