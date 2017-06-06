## input is SR 

sr = ## matrix of first column rescuer and second column vulnerable gene  
setwd("/cbcb/project2-scratch/jooslee/srescues/pancancer/cox.du/")
load("/cbcb/project2-scratch/jooslee/prob.new/prob.TCGA.mutation.correct.RData")
pancancer = prob;
rm(prob)
numGenes = nrow(pancancer$scna)
numSamples = ncol(pancancer$scna)
source("~/shortcuts/srescues/find.du.new.R")
# pval.scna 3:11 is for depletion and 12:20 for enrichment
t1 = apply(pval.scna[,c(3,14) ] ,1 , min)
t2 = apply(pval.mRNA[,c(3,14) ] ,1 , min)
molScore = ifelse(t1< 0.05 & t2 < 0.05 , t1*t2, 1 )
molecular.screen.fdr = p.adjust(molScore, method="fdr")
mol.screen.thr = ### 
molecular.screen = molScore < mol.screen.thr 


clinical.beta = cox.du.mRNA[,c(3,9)]
clinical.beta[,2] = - clinical.beta[,2] ### beta negative are desirable from here on 
clinical.p = cox.du.mRNA[,c(7,13)]
class(clinical.p) = "numeric"

clinical.p.mRNA = ifelse(clinical.beta[,1] < 0 & clinical.beta[,2] < 0, clinical.p[,1]* clinical.p[,2], 1)
clinical.fdr.mRNA = p.adjust(clinical.p.mRNA, method="fdr")


clinical.beta = cox.du.scna[,c(3,9)]
clinical.beta[,2] = - clinical.beta[,2] ### beta negative are desirable from here on 
clinical.p = cox.du.scna[,c(7,13)]
class(clinical.p) = "numeric"

clinical.p.scna = ifelse(clinical.beta[,1] < 0 & clinical.beta[,2] < 0, clinical.p[,1]* clinical.p[,2], 1)
clinical.fdr.scna = p.adjust(clinical.p.scna, method="fdr")

clincial.fdr.thr = ###

clinical.screen = which( clinical.fdr.mRNA <clincial.fdr.thr &  clinical.fdr.scna<clincial.fdr.thr ) 



### shRNA screen 
shrna.fdr.thr = ##
### Kuyuon you want the rescuer to be causal gene  
shRNA1.screen = rowSums(shrna.wil.mat[,c(2,4)] < pthr) >= 2
shRNA2.screen = rowSums(shrna.wil.bc.mat[,c(2,4)] < pthr) >=2
shRNA.screen = shRNA1.screen | shRNA2.screen 

phylo.screen.thr = ###usually we use .1 
phylogenetic.screen =( phylo.score <quantile(as.numeric(pp.rand),phylo.screen.thr,na.rm=T))  

###Kuyuon : phylogenetic screen is usually very punishing. If you don't get interactions while applying it. You can remove it and check. 
#

