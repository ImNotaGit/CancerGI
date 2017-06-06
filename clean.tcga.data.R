source("/cbcb/project2-scratch/kycheng/GI/source.goldenset.R")
load("/cbcb/project2-scratch/jooslee/prob.new/prob.TCGA.extended.qnorm.RData")
pancancer <- prob
rm(prob)
age = ifelse(is.na(pancancer$age), mean(pancancer$age,na.rm=T), pancancer$age)
race = ifelse(is.na(pancancer$race),"unknown", pancancer$race)
gii =  calculate.genomic.instability(pancancer$scna)
surv.all = data.table(pancancer$survival, types= pancancer$types, age = qnorm.array(age), sex = pancancer$sex, race = race, gii = gii)
setnames(surv.all, 1:2, c("time","status"))
surv.all$status = ifelse(surv.all$status==1,0,1)
pancancer$surv.strata = surv.all
save(pancancer, file="tcga.RData")

pancancer.surv.strata <- pancancer$surv.strata
save(pancancer.surv.strata, file="tcga.surv.strata.RData")
genes <- pancancer$genes
save(genes, file="tcga.genes.RData")
patients <- tolower(pancancer$samples)
save(patients, file="tcga.patients.RData")
pancancer.normq2 <- pancancer[c("scna.norm","mRNA.norm","scnaq2","mRNAq2")]
save(pancancer.normq2, file="tcga.normq2.RData")
scnaq2 <- pancancer$scnaq2
mRNAq2 <- pancancer$mRNAq2
save(scnaq2, mRNAq2, file="tcga.q2.RData")
scna.norm <- pancancer$scna.norm
mRNA.norm <- pancancer$mRNA.norm
save(scna.norm, mRNA.norm, file="tcga.norm.RData")
