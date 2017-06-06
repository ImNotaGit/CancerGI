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
allsr <- matrix(c(rep(1:n.genes, each=n.genes), rep(1:n.genes, n.genes)), ncol=2)

block.size <- nrow(allsr) %/% 4
block.size.last <- nrow(allsr) %% 4 + block.size
for (i in 1:4) {
  if (i!=4) {
    bs <- block.size
  } else bs <- block.size.last
  inds <- seq((i-1)*block.size+1, length.out=bs)
  scna <- hypergeometricTestPair(scnaq=scnaq2, pairs=allsr[inds,])
  save(scna, file=paste0("scna",i,".RData"))
  rm(scna)
  gc()
  mrna <- hypergeometricTestPair(scnaq=mRNAq2, pairs=allsr[inds,])
  save(mrna, file=paste0("mrna",i,".RData"))
  rm(mrna)
  gc()
}
