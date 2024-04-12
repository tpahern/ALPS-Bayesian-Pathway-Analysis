norm <- read.table('../../results/norm/testing-sampler-tree-parameters.txt',sep="\t",header=T)
> summary(norm$distance)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
0.03846 0.20000 0.33330 0.45240 0.50000 1.00000
> maxpsi <- 3
> npsi <- 100 # number of psi's to explore
> psi=(maxpsi*1:npsi)/npsi
>
> cpsi<-array(dim=npsi)
> for (i in 1:npsi) {
+       cpsi[i] <- sum(exp(-psi[i]*norm$distance))/nrow(norm)
+ }
> cpsi