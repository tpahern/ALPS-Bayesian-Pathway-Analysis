# Simulate complex pathway for ALPS example data

library(ggplot2)
library(survival)

wd <- "C:\\Users\\tomah\\OneDrive\\Documents\\GitHub\\ALPS-Bayesian-Pathway-Analysis\\simpath"
setwd(wd)

# import observed (but scrambled) ProBe CaRe cohort genotype data
simgeno <- read.csv("cohort_sim1.csv", header=TRUE)

attach(simgeno)  


# simpath.png in working directory shows the tamoxifen paths we want to simulate
  # genes & variants involved: 
    # CYP2C9 -  cyp7
    # CYP2C19 - cyp4
    # ABCG2 - abc9
    # SULT1A1 - sul3

# Will need six theta_lambda parameters and one beta parameter to govern the pathway
# (cyp7, cyp4) gets two thetas
# (abc9, sul3) gets two thetas
# ([cyp7,cyp4],[abc9,sul3]) gets two thetas

# specify thetas

  # (cyp7,cyp4)
  theta1 <- 0
  theta2 <- 0

  # (abc9,sul3)
  theta3 <- 1
  theta4 <- 1

  # (cyp7,cyp4)(abc9,sul3)
  theta5 <- 0.5
  theta6 <- 0.5

# create nested pathway structures
z1 <- cyp7*theta1 + cyp4*theta2 + (1-theta1-theta2)*cyp7*cyp4
z2 <- abc9*theta3 + sul3*theta4 + (1-theta3-theta4)*abc9*sul3
simgeno$z <- z1*theta5 + z2*theta6 + (1-theta5-theta6)*z1*z2

# simulate survival data based on z, beta
  # simulation parameters
  beta <- -3 # the overall pathway coefficient
  lambda.t <- 0.01  #baseline hazard
  lambda.c <- 0.005 #censoring hazard

linpred <- exp(beta*simgeno$z)
scale.t <- linpred*lambda.t

  set.seed(557)
  t <- rweibull(n=nrow(simgeno), shape=1, scale=scale.t) #time to event
  c <- rweibull(n=nrow(simgeno), shape=1, lambda.c) #time to censoring
  simgeno$time <- round(100000*(pmin(t,c)),2) #follow-up time (whichever comes first)
  
  #status variable
  censored <- ifelse(c < t,1,0)
  simgeno$event <- 1-censored

# check distribution of person-time
ggplot(simgeno, aes(x=time)) + geom_histogram()

mdl_marginal_main <- coxph(Surv(time,event) ~ 
                             (abc1	+ abc2 + abc3	+ abc4 + abc5	+ abc6 + abc7	+
                                abc8	+ abc9 + abc10 + cyp1	+ cyp2 + cyp3	+ cyp4 +
                                cyp5	+ cyp6 + cyp7	+ cyp8 + cyp9	+ cyp10 + cyp11 +
                                cyp12+ cyp13 + cyp14 + sul1 + sul2	+ sul3 + sul4	+
                                sul5	+ ugt1 + ugt2	+ ugt3), 
                           data=simgeno)
summary(mdl_marginal_main)

# write out simulated pathway data for use with ALPS
write.csv(simgeno,file="simpath_cohort_scram.csv")