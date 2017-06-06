library(parallel)
library(doMC)
library(foreach)
ncores=23L
registerDoMC(ncores)

load("/cbcb/project2-scratch/kycheng/GI/tcga.genes.RData") # genes

## simplify the data for shRNA screening: run once then comment out as record
# dataset 1
#load("prob_essentiality.Dan.RData") # prob
#ess <- list(genes=prob$genes, scna=prob$scna, mRNA=prob$mRNA, ess=prob$ess)
#save(ess, file="my.dan.RData")
# dataset 2
#load("sayad.RData") # dat
#ess <- list(genes=dat$genes, scna=dat$scna, mRNA=dat$mRNA, ess=dat$mat)
#save(ess, file="my.sayad.RData")
# I include some datasets for the screens in another version of the pipeline too:
# old screen 1
#load("achilles.old.RData") # prob # Cheung et al. PNAS (2011).
#ess <- list(genes=prob$genes, scna=prob$scna, mRNA=prob$mRNA, ess=prob$mat)
#save(ess, file="my.achilles.old.RData")
# old screen 2
#load("achilles.new.RData") # prob # Cowley et al. Sci. Data. (2014).
#ess <- list(genes=prob$genes, scna=prob$scna, mRNA=prob$mRNA, ess=prob$mat)
#save(ess, file="my.achilles.new.RData")
# old screen 3
#load("marcotte.old.RData") # prob # Marcotte et al. Cancer Discov (2012).
#ess <- list(genes=prob$genes, scna=prob$scna, mRNA=prob$mRNA, ess=prob$mat)
#save(ess, file="my.marcotte.old.RData")
# old screen 4
#load("marcotte.new.RData") # prob # Marcotte et al. Cell (2016).
#ess <- list(genes=prob$genes, scna=prob$scna, mRNA=prob$mRNA, ess=prob$mat)
#save(ess, file="my.marcotte.new.RData")

ess.dat1 <- "/cbcb/project2-scratch/kycheng/GI/my.dan.RData"
ess.dat2 <- "/cbcb/project2-scratch/kycheng/GI/my.sayad.RData"
ess.dat3 <- "/cbcb/project2-scratch/kycheng/GI/my.achilles.old.RData"
ess.dat4 <- "/cbcb/project2-scratch/kycheng/GI/my.achilles.new.RData"
ess.dat5 <- "/cbcb/project2-scratch/kycheng/GI/my.marcotte.old.RData"
ess.dat6 <- "/cbcb/project2-scratch/kycheng/GI/my.marcotte.new.RData"

wilcox.test.na = function(x,y, alternative1, paired=FALSE) {
    tryCatch(
        wilcox.test(x, y, alternative=alternative1, paired=paired)$p.value,
        error = function(e) NA
        )
}

pair.ess.wil = function(pair, ess=ess)
{
    out = rep(NA, 8)
    if(sum(is.na(pair)) == 0){
      aa = ess$mRNA[pair[1],]
      out[1]= wilcox.test.na(ess$ess[pair[2],aa <= median(aa,na.rm=T)],ess$ess[pair[2],aa >median(aa,na.rm=T)],alternative1="less" )
      aa = ess$mRNA[pair[2],]
      out[2]= wilcox.test.na(ess$ess[pair[1],aa <= median(aa,na.rm=T)],ess$ess[pair[1],aa >median(aa,na.rm=T)],alternative1="less" )
      aa = ess$scna[pair[1],]
      out[3]= wilcox.test.na(ess$ess[pair[2],aa <= median(aa,na.rm=T)],ess$ess[pair[2],aa >median(aa,na.rm=T)],alternative1="less" )
      aa = ess$scna[pair[2],]
      out[4]= wilcox.test.na(ess$ess[pair[1],aa <= median(aa,na.rm=T)],ess$ess[pair[1],aa >median(aa,na.rm=T)],alternative1="less" )

      aa = ess$mRNA[pair[1],]
      out[5]= wilcox.test.na(ess$ess[pair[2],aa <= median(aa,na.rm=T)],ess$ess[pair[2],aa >median(aa,na.rm=T)],alternative1="greater" )
      aa = ess$mRNA[pair[2],]
      out[6]= wilcox.test.na(ess$ess[pair[1],aa <= median(aa,na.rm=T)],ess$ess[pair[1],aa >median(aa,na.rm=T)],alternative1="greater" )
      aa = ess$scna[pair[1],]
      out[7]= wilcox.test.na(ess$ess[pair[2],aa <= median(aa,na.rm=T)],ess$ess[pair[2],aa >median(aa,na.rm=T)],alternative1="greater" )
      aa = ess$scna[pair[2],]
      out[8]= wilcox.test.na(ess$ess[pair[1],aa <= median(aa,na.rm=T)],ess$ess[pair[1],aa >median(aa,na.rm=T)],alternative1="greater" )
  }
  out
}

calc.rawp.ess <- function(gpis) {
  lapply(paste0("ess.dat", 1:6), function(dat) {
    load(get(dat)) # ess
    sr.ess <- cbind(match(genes[gpis[,1]], ess$genes), match(genes[gpis[,2]], ess$genes))
    res <- foreach(x=1:nrow(sr.ess), .inorder=TRUE, .combine=rbind) %dopar% pair.ess.wil(sr.ess[x,], ess=ess)
    res <- as.data.table(res)
    res[, c("rn","vn","ri","vi"):=list(genes[gpis[,1]], genes[gpis[,2]], gpis[,1], gpis[,2])]
    setcolorder(res, c("rn","vn","ri","vi",names(res)[1:(length(res)-4)]))
    res
  })
}
