#build DU SRs
if(!exists("preprocessed")) preprocessed =FALSE
if(!exists("dont.run.other.screen")) dont.run.other.screen =FALSE
if(!preprocessed){
#build DU SRs
source("~/shortcuts/srescues/source.goldenset.R")
library(parallel)
numGenes = length(pancancer$genes)
numSamples = length(pancancer$sample)
mRNA.norm = mclapply(1:numSamples, function(tt) qnorm.array(pancancer$mRNA[,tt]),mc.cores=64)
scna.norm = mclapply(1:numSamples, function(tt) qnorm.array(pancancer$scna[,tt]),mc.cores=64)
mRNA.norm = t(do.call(rbind, mRNA.norm))
scna.norm = t(do.call(rbind, scna.norm))

mRNA.norm = mclapply(1:numGenes, function(tt) qnorm.array(mRNA.norm[tt,]),mc.cores=64)
scna.norm = mclapply(1:numGenes, function(tt) qnorm.array(scna.norm[tt,]),mc.cores=64)
mRNA.norm = do.call(rbind, mRNA.norm)
scna.norm = do.call(rbind, scna.norm)
scnaq2 = pancancer$scnaq2
mRNAq2 = pancancer$mRNAq2


pancancer$mRNA.norm = mRNA.norm 
pancancer$scna.norm = scna.norm
pancancer$age = ifelse(is.na(pancancer$age), mean(pancancer$age,na.rm=T), pancancer$age)
pancancer$race = ifelse(is.na(pancancer$race),"unknown", pancancer$race)
pancancer$sex = ifelse(is.na(pancancer$sex),"FEMALE", pancancer$sex)
gii =  calculate.genomic.instability(pancancer$scna)
surv.all = data.table(pancancer$survival, types= pancancer$types, age = qnorm.array(pancancer$age), sex = pancancer$sex, race = pancancer$race, gii = gii)
setnames(surv.all, 1:2, c("time","status"))
surv.all$status = ifelse(surv.all$status==1,0,1)
pancancer$surv.strata = surv.all
pancancer$types = pancancer$type 
}
### clinical 
f2.list = c(  "strata(sex)",  "age", "strata(race)", "gii")
f1.list = "strata(types)"
cox.mRNA = mclapply(1:nrow(sr), function(tt) cox.pair.du.strata.controlled(sr[tt,], prob = pancancer, use.mRNA=T, f2.list = f2.list, f1.list = f1.list),mc.cores=64)
cox.du.mRNA = do.call(rbind, cox.mRNA)
cox.du.mRNA=cbind(sr,cox.du.mRNA)

cox.scna = mclapply(1:nrow(sr), function(tt) cox.pair.du.strata.controlled(sr[tt,], prob = pancancer, use.mRNA=F, f2.list = f2.list, f1.list = f1.list),mc.cores=64)
# x.scna = mclapply(1:128, function(tt) cox.pair.du.strata.controlled(sr[tt,], prob = pancancer, use.mRNA=F, f2.list = f2.list),mc.cores=64)
cox.du.scna = do.call(rbind, cox.scna)
cox.du.scna=cbind(sr,cox.du.scna)

## molecular 
library(Rcpp)
sourceCpp("~/shortcuts/srescues/HyperGeometricTest.pair.cpp", rebuild=T)
pval.scna1 = hypergeometricTestPair(scnaq= pancancer$scnaq2, pairs=sr)
pval.scna.up = hypergeometricTestPair(scnaq= pancancer$scnaq2, pairs=sr,lowerTail=0)
pval.scna = cbind(sr[,1:2], pval.scna1, pval.scna.up)

pval.mRNA1 = hypergeometricTestPair(scnaq= pancancer$mRNAq2, pairs=sr)
pval.mRNA.up = hypergeometricTestPair(scnaq= pancancer$mRNAq2, pairs=sr,lowerTail=0)
pval.mRNA = cbind(sr[,1:2], pval.mRNA1, pval.mRNA.up)

##shrna 
#### calculate correlation with essentiality 
if(!dont.run.other.screen){
load("/cbcb/project2-scratch/jooslee/srescues/shrna/prob_essentiality.Dan.RData")
ess=prob
rm(prob)
sr.ess = cbind( match(pancancer$genes[sr[,1]],ess$genes), match(pancancer$genes[sr[,2]],ess$genes))
require(doMC)
require(foreach)
registerDoMC(cores = 64)
outbed = foreach(x = 1:nrow(sr.ess), .inorder=T) %dopar% pair.ess.wil(sr.ess[x,], ess = ess)
shrna.wil.mat = do.call(rbind, outbed)


ess.bc = local({load("/cbcb/project2-scratch/jooslee/prob.new/sayad.RData"); environment()})
ess.bc=ess.bc$dat
ess.bc$ess = ess.bc$mat
sr.ess.bc = cbind( match(pancancer$genes[sr[,1]],ess.bc$genes),
	match(pancancer$genes[sr[,2]],ess.bc$genes))
outbed = foreach(x = 1:nrow(sr.ess.bc), .inorder=T) %dopar% pair.ess.wil(sr.ess.bc[x,],ess=ess.bc)
shrna.wil.bc.mat = do.call(rbind, outbed)

## phylogenetic
sr.gene.all = data.table(rescuer = pancancer$genes[sr[,1]], vulnerable = pancancer$genes[sr[,2]])
phylo.score = phylo.profile(sr.gene.all )

sr.rand=data.table(rescuer=pancancer$genes[sample(numGenes,100000,replace=T)],vulnerable=pancancer$genes[sample(numGenes,100000,replace=T)])
na.inx=which(sr.rand$rescuer!=sr.rand$vulnerable)
sr.rand=sr.rand[na.inx,]
pp.rand=phylo.profile(sr.rand)
}
