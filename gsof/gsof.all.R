#!/opt/local/stow/R-3.3.1/bin/Rscript

library(data.table)
library(Rcpp)
library(parallel)
library(doMC)
library(foreach)
ncores <- 8L
registerDoMC(ncores)

setwd("/cbcb/project2-scratch/kycheng/cancer_SR_drug/gsof")
load("../tcga.genes.RData")
load("../tcga.q2.RData")
sourceCpp("../HyperGeometricTest.pair.cpp")


### ----gSoF----

n.genes <- length(genes)
block.size <- n.genes %/% 4
block.size.last <- n.genes %% 4 + block.size

#mcopt <- list(preschedule=FALSE)
for (i in 1:4) {
  if (i!=4) {
    bs <- block.size
  } else bs <- block.size.last
  xs <- seq((i-1)*block.size+1, length.out=bs)
  
  scna <- foreach(x=xs, .inorder=TRUE, .combine=rbind) %dopar% {#, .options.multicore=mcopt) %dopar% {
    hypergeometricTest(scnaq2[x,], scnaq2)## conduct hypergeometric test for all possible 9 functional states of a gene pair.
  }
  save(scna, file=paste0("scna",i,".RData"))
  rm(scna)
  gc()
  
  mrna <- foreach(x=xs, .inorder=TRUE, .combine=rbind) %dopar% {#, .options.multicore=mcopt) %dopar% {
    hypergeometricTest(mRNAq2[x,], mRNAq2)## conduct hypergeometric test for all possible 9 functional states of a gene pair.
  }
  save(mrna, file=paste0("mrna",i,".RData"))
  rm(mrna)
  gc()
}
