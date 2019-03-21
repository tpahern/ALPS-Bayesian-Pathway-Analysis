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
source(file=paste0(wd,"alps2-test.R"))

# Read in Dataset
ds <- read.csv(paste0(wd,'survsim.csv'),header=T)

# Read in SNP-level prior specification
all <- read.table(file=paste0(wd,"prior-forest-snps.txt"),sep="\t",header=T)
prior.forest <- subset(all,select=c("snp1","snp2"));

# Config dataset for ALPS
time.var <- ds$time
status.var <- ds$event

# Note that it is called dos for genotype dosages. 
dos <- ds;
dos$labid <- NULL
dos$time <- NULL
dos$event <- NULL

dos <- t(as.matrix(dos))

#-------------------------------------------------------------
# Initialize ALPS and run
#-------------------------------------------------------------
load(file=paste0(wd,'cpsi.RData'))

# start at a random spot in the prior forest
spot <- prior.forest[sample(nrow(prior.forest),1),]
spot.tree <- paste0("(",spot$snp1,",",spot$snp2,");")
curtree <- read.tree(text=spot.tree)

#curtree <- read.tree(text="9((gdosecyp8,gdoseabc9),(gdosesul1,gdosecyp6));")


# Tom, change iterations.  I would see how long 1000 takes, then scale up. Goal is to get 100,000.
fitalps(iter=10,initipsi=19,curtree,normpotts=FALSE,prioronly=FALSE)





