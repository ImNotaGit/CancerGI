##### ISLE source code #####
##### Joo Sang Lee, Avinash Das, Eytan Ruppin #####

############# load libraries
library(data.table)
library(Rcpp)
library(parallel)
library(survival)

require(doMC)
require(foreach)
registerDoMC(cores = 64)

####### load TCGA data ##########
load("TCGA.RData")

####### The mRNA expression (RNAseqV2) and patients' clinical characteristics were downloaded from TCGA Data Portal,
### and copy number data (Gistic values) was downloaded from Broad Institute's Firehose (https://gdac.broadinstitute.org/)
### on June 27, 2015. The TCGA data includes the following data fields:
### genes         19001 -none-     character: protein coding genes' symbols
### entrez        19001 -none-     numeric  : protein coding genes' entrez ID
### samples        8749 -none-     character: TCGA sample barcodes
### types          8749 -none-     character: cancer types 
### mRNA      166239749 -none-     numeric  : a matrix of gene expression (RNAseq) with genes in the row and samples in the column
### scna      166239749 -none-     numeric  : a matrix of SCNA (Gistic values) with genes in the row and samples in the column
### mRNAq2    166239749 -none-     numeric  : a matrix of mRNA expression equally divided into 3 groups (low:0, middle:1, high:2) using 1/3, 2/3-quantile across the samples in each cancer type
### scnaq2    166239749 -none-     numeric  : a matrix of SCNA equally divided into 3 groups (low:0, middle:1, high:2) using 1/3, 2/3-quantile across the samples in each cancer type
### sex            8749 -none-     character: sex of patients
### age            8749 -none-     numeric  : age of patients 
### race           8749 -none-     character: ethnicity of patients
### survival      17498 -none-     numeric  : two column matrix with first column as survival time (in days) and second column as cencering (death=0, alive=1)
### stage          8749 -none-     character: stage of the tumor 
### mRNA.norm 166239749 -none-     numeric  : quantile-normalized of mRNA expression values
### scna.norm 166239749 -none-     numeric  : quantile-normalized of SCNA values
### surv.dt           2 data.table list		: survival matrix for Cox regression (first column: survival time (in days), second column: censoring (death=1, alive=0))


############# function for quantile-normalization
qnorm.array <- function(mat)     
{
	mat.back = mat 
	mat = mat[!is.na(mat)]
    mat = rank(mat, ties.method = "average");
    mat = qnorm(mat / (length(mat)+1));
    mat.back[!is.na(mat.back)] = mat 
    mat.back
}

############# store variables for following analysis
genes=prob$genes
numGenes=length(genes)
numSamples=length(prob$samples)
cancerType = prob$types;typeInx = as.numeric(as.factor(cancerType))
typeNum= length(unique(typeInx ));levels(typeInx) = seq(typeNum)-1
mRNA.norm = prob$mRNA.norm
scna.norm = prob$scna.norm
surv.dt = prob$surv.dt
surv.type=as.numeric(as.factor(prob$types))
age=qnorm.array(prob$age)
race=prob$race
sex=prob$sex

##########################################################################################
################### Collecting P-values from the 4 ISLE screenings #######################
##########################################################################################

##########################################################################################
###### Step 1. gSoF screening ############################################################
##########################################################################################
sourceCpp("binomialTest.cpp")
numGenes = nrow(scna)
numSamples = ncol(scna)
scnaq=prob$scnaq2
mRNAq=prob$mRNAq2
xs=seq(numGenes)

########## molecular scrneeing (binomial test) for mRNA:
########## outbed[,1:2]: gene indices, outbed[,3]: the significance of gene A-gene B low depletion
########## negative p-values -> depletion, positive p-values -> enrichment
outbed = foreach(x = xs, .inorder=T, .combine=rbind) %dopar%{
	mRNAcurr = mRNAq[x,] 
	mRNAcurr[is.na(mRNAcurr)] = 1
	out = binomialTest(mRNAcurr,mRNAq)
	out2=cbind(x,seq(numGenes),out)
	print(x)
	return(out2)
}
inx.mRNA=which(abs(outbed[,3])<0.1 & outbed[,3]<0 )
sl.mRNA=outbed[,1:3]
sl.id=outbed[,1:2] # store gene indices

########## molecular scrneeing (binomial test) for SCNA:
########## outbed[,1:2]: gene indices, outbed[,3]: the significance of gene A-gene B low depletion
########## negative p-values -> depletion, positive p-values -> enrichment
outbed = foreach(x = xs, .inorder=T, .combine=rbind) %dopar%{
	scnacurr = scnaq[x,] 
	scnacurr[is.na(scnacurr)] =1
	out = binomialTest(scnacurr,scnaq)
	out2=cbind(x,seq(numGenes),out)
	print(x)
	return(out2)
}
inx.scna=which(abs(outbed[,3])<0.1 & outbed[,3]<0 )
sl.scna=outbed[,1:3]

########## obtain thresholds for empirical false discovery correction
q.mRNA=quantile(abs(sl.mRNA[,3]),0.1) # empirical threshold for mRNA screening (10% of FDR correction)
q.scna=quantile(abs(sl.scna[,3]),0.1) # empirical threshold for scna screening (10% of FDR correction)

########## transform negative p-values (for depletion) to be positive, and positive p-values (enrichment) to be 1
sl.mRNA[,3]=ifelse(sl.mRNA[,3]<0,-sl.mRNA[,3],1)
sl.scna[,3]=ifelse(sl.scna[,3]<0,-sl.scna[,3],1)

########## Select pairs that pass Step 1: gSoF screening
i.mol=sl.mRNA[,3]<q.mRNA & sl.scna[,3]<q.scna
sl.mol=cbind(sl.id[i.mol,],sl.mRNA[i.mol,3],sl.scna[i.mol,3])
sl.id=sl.id[i.mol,]

##########################################################################################
###### Step 2. Clinical screening ########################################################
##########################################################################################

########## function for Cox regression 
cox.pair.sl = function(pair,use.mRNA=F)
{
	if(use.mRNA){
		g1 = mRNA.norm[pair[1],]  
		g2 = mRNA.norm[pair[2],]
		f1 = prob$mRNAq2[pair[1],]  
		f2 = prob$mRNAq2[pair[2],]  
	}else{
		g1 = scna.norm[pair[1],]  
		g2 = scna.norm[pair[2],]  
		f1 = prob$scnaq2[pair[1],]  
		f2 = prob$scnaq2[pair[2],]  
	}
	cov = ifelse(f1 == 0 & f2 ==0,1,0 )
	dt1 = data.frame(cbind(surv.dt, surv.type, cov, cbind(g1 , g2)))
	cntrl.out = coxph(Surv(time,status) ~ g1 + g2 + strata(surv.type), data=dt1)
	ll1 = cntrl.out$loglik[2]
	cox.out = coxph(Surv(time,status) ~ cov + g1 + g2 + strata(surv.type), data=dt1)
	ll2 = cox.out$loglik[2]
	aa  = summary(cox.out)
	uncntrl.out = coxph(Surv(time,status) ~ cov + strata(surv.type), data=dt1)	
	bb  = summary(uncntrl.out)	
	c(bb$coefficients["cov",], aa$coefficients["cov",], ll2 - ll1)
} 
########## function for Cox regression with controlling for age, gender, race
cox.pair.sl2 = function(pair,use.mRNA=F)
{
	if(use.mRNA){
		g1 = mRNA.norm[pair[1],]  
		g2 = mRNA.norm[pair[2],]
		f1 = prob$mRNAq2[pair[1],]  
		f2 = prob$mRNAq2[pair[2],]  
	}else{
		g1 = scna.norm[pair[1],]  
		g2 = scna.norm[pair[2],]  
		f1 = prob$scnaq2[pair[1],]  
		f2 = prob$scnaq2[pair[2],]  
	}
	cov = ifelse(f1 == 0 & f2 ==0,1,0 )
	dt1 = data.frame(cbind(surv.dt, surv.type, cov, cbind(g1 , g2), age, sex, race))
	cntrl.out = coxph(Surv(time,status) ~ g1 + g2 + age + strata(surv.type, sex, race), data=dt1)
	ll1 = cntrl.out$loglik[2]
	cox.out = coxph(Surv(time,status) ~ cov + g1 + g2 + age + strata(surv.type, sex, race), data=dt1)
	ll2 = cox.out$loglik[2]
	aa  = summary(cox.out)
	uncntrl.out = coxph(Surv(time,status) ~ cov + age + strata(surv.type,sex,race), data=dt1)	
	bb  = summary(uncntrl.out)	
	c(bb$coefficients["cov",], aa$coefficients["cov",], ll2 - ll1)
}

#### clinical screening (based on mRNA)
	cox.mRNA = mclapply(1:nrow(sl.id), function(tt) cox.pair.sl(sl.id[tt,], use.mRNA=T),mc.cores=64)
	cox.sl.mRNA = t(do.call(cbind, cox.mRNA))
cox.sl.mRNA=cbind(sl,cox.sl.mRNA)

#### clinical screening (based on SCNA)
	cox.scna = mclapply(1:nrow(sl.id), function(tt) cox.pair.sl(sl.id[tt,], use.mRNA=F),mc.cores=64)
	cox.sl.scna = t(do.call(cbind, cox.scna))
cox.sl.scna=cbind(sl,cox.sl.scna)

#### clinical screening to control for age, sex, race (based on mRNA)
	cox.mRNA = mclapply(1:nrow(sl.id), function(tt) cox.pair.sl2(sl.id[tt,], use.mRNA=T),mc.cores=64)
	cox.sl.mRNA2 = t(do.call(cbind, cox.mRNA))
cox.sl.mRNA2=cbind(sl,cox.sl.mRNA2)

#### clinical screening to control for age, sex, race (based on SCNA)
	cox.scna = mclapply(1:nrow(sl.id), function(tt) cox.pair.sl2(sl.id[tt,], use.mRNA=F),mc.cores=64)
	cox.sl.scna2 = t(do.call(cbind, cox.scna))
cox.sl.scna2=cbind(sl,cox.sl.scna2)

########## Select pairs that pass Step 2: Clinical screening
pcox=p.adjust(c(cox.sl.mRNA[,7],cox.sl.mRNA[,12]),"BH")
i.cox1=cox.sl.mRNA[,3]<0 & pcox[1:nrow(cox.sl.mRNA)]<0.1 & cox.sl.mRNA[,8]<0 & pcox[1:nrow(cox.sl.mRNA)+nrow(cox.sl.mRNA)]<0.1

pcox=p.adjust(c(cox.sl.scna[,7],cox.sl.scna[,12]),"BH")
i.cox2=cox.sl.scna[,3]<0 & pcox[1:nrow(cox.sl.scna)]<0.1 & cox.sl.scna[,8]<0 & pcox[1:nrow(cox.sl.scna)+nrow(cox.sl.scna)]<0.1

pcox=p.adjust(c(cox.sl.mRNA2[,7],cox.sl.mRNA2[,12]),"BH")
i.cox3=cox.sl.mRNA2[,3]<0 & pcox[1:nrow(cox.sl.mRNA2)]<0.1 & cox.sl.mRNA2[,8]<0 & pcox[1:nrow(cox.sl.mRNA2)+nrow(cox.sl.mRNA2)]<0.1

pcox=p.adjust(c(cox.sl.scna2[,7],cox.sl.scna2[,12]),"BH")
i.cox4=cox.sl.scna2[,3]<0 & pcox[1:nrow(cox.sl.scna2)]<0.1 & cox.sl.scna2[,8]<0 & pcox[1:nrow(cox.sl.scna2)+nrow(cox.sl.scna2)]<0.1

i.cox=which((i.cox1 | i.cox2) & (i.cox3 | i.cox4) )
sl.id=sl.id[i.cox,]
sl.cox=sl.cox[i.cox,]

##########################################################################################
###### Step 3. Phenotypic screening ######################################################
##########################################################################################
### We collected the following 4 large-scale shRNA essentiality screening datasets for step 3, 
### and run the screening for all possible gene pairs
load("achilles.old.RData") # Cheung et al. PNAS (2011).
load("achilles.new.RData") # Cowley et al. Sci. Data. (2014).
load("marcotte.old.RData") # Marcotte et al. Cancer Discov (2012).
load("marcotte.new.RData") # Marcotte et al. Cell (2016).

### this table shows the data fields included in the first dataset (Cheung et al.) as an example,
### the remaining 3 datasets have a similar data structure.
###          Length  Class      Mode     
###genes       19001 -none-     character : gene symbols of protein coding genes
###entrez      19001 -none-     numeric   : entrez IDs of protein coding genes
###samples       102 -none-     character : sample IDs (as coded in the dataset)
###celllines     102 -none-     character : names of cell lines
###types         102 -none-     character : cancer types of cell lines
###mRNA      1938102 -none-     numeric   : matrix of mRNA expression with genes in rows, cell lines in columns
###mRNAq     1938102 -none-     numeric   : matrix of mRNA expression equally divided into 3 groups (low:0, middle:1, high:2) using 1/3, 2/3-quantile across cell lines
###scna      1938102 -none-     numeric   : matrix of SCNA expression with genes in rows, cell lines in columns
###scnaq     1938102 -none-     numeric   : matrix of SCNA expression equally divided into 3 groups (low:0, middle:1, high:2) using 1/3, 2/3-quantile across cell lines
###mat           102 data.frame list      : matrix of shRNA essentiality scores with genes in rows, cell lines in columns (the lower the more essential)

### function to remove NA values in essentiality score matrix
na.solve=function(mat){ 
	medv=apply(mat,1,median,na.rm=T)
	medv[is.na(medv)]=median(medv,na.rm=T)
	mat=do.call(rbind,lapply(seq(nrow(mat)),function(i) {v=mat[i,];v[is.na(v)]=medv[i];v}))
	mat
}

### map gene indices, here the screening is asymmetric because we are using different types of data for each gene - essentiality and expression/SCNA
### to consider this asymmetry, we map the indices so that the screening for A-B and B-A is stored in the same row
### this makes the indicing compatible with the previous 2 screenings
sl.id=cbind(rep(seq(numGenes),each=numGenes),rep(seq(numGenes),numGenes))
ix=which(sl.id[,1]<sl.id[,2])
sl.id1=sl.id[ix,]
usr=paste(sl.id[,1],sl.id[,2])
usr1=paste(sl.id1[,1],sl.id1[,2])
usr2=paste(sl.id1[,2],sl.id1[,1])
idz1=match(usr1,usr)
idz2=match(usr2,usr)
sourceCpp("essentialityTest.cpp")

### the following code shows how to obtain the phenotypic significance of all gene pairs
### for the first dataset (Cheung et al.) as an example
#########################################################################################################
### one can repeat the process for the remaining 3 datasets to obtain p.ach.new, p.mar.old, p.mar.new ###
#########################################################################################################

prob$mat=do.call(cbind,lapply(seq(ncol(prob$mat)),function(i) as.numeric(as.matrix(prob$mat[,i]))))
prob$mat=na.solve(prob$mat)
nSmp=length(prob$samples)
scnatest=prob$mRNA
	scna1=scnatest> apply(scnatest, 1, quantile, probs = 1/2,  na.rm = TRUE)
	prob$mRNAq3=scna1*2
scnatest=prob$scna
	scna1=scnatest> apply(scnatest, 1, quantile, probs = 1/2,  na.rm = TRUE)
	prob$scnaq3=scna1*2
mRNAq3=prob$mRNAq3
scnaq3=prob$scnaq3
mRNAq3[is.na(mRNAq3)]=1
scnaq3[is.na(scnaq3)]=1
mat=prob$mat
p.ach.old=essentialityTestAll(t(scnaq3), t(mRNAq3), t(mat), 1)
p.ach.old=cbind(sl.id1,p.ach.old[idz1,],p.ach.old[idz2,])
p.ach.old[p.ach.old==-1000]=NA

# p.ach.old[,3]: when gene 1 is more active, gene 2 has higher essentiality score (less essential) - SCNA
# p.ach.old[,4]: when gene 1 is more active, gene 2 has higher essentiality score (less essential) - mRNA
# p.ach.old[,5]: when gene 2 is more active, gene 1 has higher essentiality score (less essential) - SCNA
# p.ach.old[,6]: when gene 2 is more active, gene 1 has higher essentiality score (less essential) - mRNA

### select the phenotypic screening results of relevant gene pairs given by sl.cox
usr=paste(sl.cox[,1],sl.cox[,2])
usr.ess=paste0(p.ach.old[,1],"&",p.ach.old[,2])
ix1=match(usr,usr.ess)
p.ach.old=p.ach.old[ix1,3:6];p.ach.new=p.ach.new[ix1,3:6]
p.mar.old=p.mar.old[ix1,3:6];p.mar.new=p.mar.new[ix1,3:6]

class(p.ach.old) <- "numeric";class(p.ach.new) <- "numeric"
class(p.mar.old) <- "numeric";class(p.mar.new) <- "numeric"

# remove p-value of 0.5 that is the case where all cell line was assigned the same essentiality (in the process of replacing NAs)
p.ach.old[p.ach.old==0.5]=NA;p.ach.new[p.ach.new==0.5]=NA
p.mar.old[p.mar.old==0.5]=NA;p.mar.new[p.mar.new==0.5]=NA

p.ess=cbind(
apply(cbind(apply(p.ach.old[,1:2],1,min,na.rm=T),apply(p.ach.old[,3:4],1,min,na.rm=T)),1,min,na.rm=T),
apply(cbind(apply(p.ach.new[,1:2],1,min,na.rm=T),apply(p.ach.new[,3:4],1,min,na.rm=T)),1,min,na.rm=T),
apply(cbind(apply(p.mar.old[,1:2],1,min,na.rm=T),apply(p.mar.old[,3:4],1,min,na.rm=T)),1,min,na.rm=T),
apply(cbind(apply(p.mar.new[,1:2],1,min,na.rm=T),apply(p.mar.new[,3:4],1,min,na.rm=T)),1,min,na.rm=T))
class(p.ess) <- "numeric"

### A-B was marked to show the phenotypic sifnificance if either mRNA or SCNA support the interaction 
### based on either conditional essentiality of A or B
p.ess=as.numeric(apply(p.ess,1,min,na.rm=T))
p.ess[is.infinite(p.ess)]=NA

########## Select pairs that pass Step 3: Phenotypic screening
p.ess=p.adjust(p.ess[,25],"BH")
i.ess=which(p.ess<0.2)
sl.ess=sl.cox[i.ess,]
sl.id=sl.ess[,1:2]

##########################################################################################
###### Step 4. Phylogenetic screening ####################################################
##########################################################################################
load("yuval.phylogenetic.profile.RData")
### The phylogenetic profile is downloaded from Yuval Tabach et al. Mol Syst Biol. (2013), Supplementary Table 1
load("feature.weight.RData")
### the feature weights are determined based on the phylogenetic tree (Ensembl database: http://useast.ensembl.org/index.html)

### function to identify phylogenetic distance of a pair of genes
phylo.profile = function(sl.gene.all){
    sl.gene1 = sl.gene.all
    sl.phylo =  cbind(match(sl.gene1[,1], phylo$genes), match(sl.gene1[,2], phylo$genes))
    featureMat = (phylo[sl.phylo[,1],-(1:3)] - phylo[sl.phylo[,2],-(1:3)])^2
    featureMat %*% t(feature.weight)
}

### calculate the phylogenetic distance of our SL pairs
pp=phylo.profile(cbind(prob$genes[sl.id[,1]],prob$genes[sl.id[,2]]))

### calculate the phylogenetic distance of random gene pairs for false discovery correction
sr.rand=cbind(prob$genes[sample(numGenes,100000,replace=T)],prob$genes[sample(numGenes,100000,replace=T)])
na.inx=which(sr.rand[,1]!=sr.rand[,2])
sr.rand=sr.rand[na.inx,]
pp.rand=phylo.profile(sr.rand)

########## Select pairs that pass Step 4: Phylogenetic screening
i.pp=which(i.mol & (sl.full[,26]<quantile(as.numeric(pp.rand),0.1,na.rm=T)))
sl.pp=sl.ess[i.pp,]

##########################################################################################
########## Final SL pairs that pass all 4 screenings #####################################
##########################################################################################

sl.final=cbind(prob$genes[sl.pp[,1]],prob$genes[sl.pp[,2]])











