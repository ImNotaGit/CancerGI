library(parallel)
library(survival)
ncores=23L

load("/cbcb/project2-scratch/kycheng/GI/tcga.genes.RData") # genes
load("/cbcb/project2-scratch/kycheng/GI/tcga.normq2.RData") # pancancer.normq2
load("/cbcb/project2-scratch/kycheng/GI/tcga.surv.strata.RData") # pancancer.surv.strata
pancancer <- pancancer.normq2
pancancer$surv.strata <- pancancer.surv.strata
rm(pancancer.normq2, pancancer.surv.strata)
gc()

qnorm.array <- function(mat)
{
  mat.back = mat 
  mat = mat.back[!is.na(mat.back)]
    mat = rank(mat,  rank, ties.method = "average");
    mat = qnorm(mat / (length(mat)+1));
    mat.back[!is.na(mat.back)] = mat 
    mat.back
}

coxph.run = function(formula.list,data1){
  fmla <- as.formula(paste(" Surv(time,status) ~ ", paste(formula.list, collapse= "+")))
  cox.out = coxph(fmla, data=data1)
  ll = cox.out$loglik[2]
  aa  = summary(cox.out)
  out = ll 
  if("cov" %in% rownames(aa$coefficients) ){
    out = c(aa$coefficients["cov", ], ll)
  }
  out
}

coxph.robust = function(data1, f1.list, f2.list=NULL){
  tryCatch( 
    coxph.run(c(f1.list, f2.list), data1),
    error = function(e)  coxph.run(f1.list, data1)
    )
}

clin.screen <- function(pair, prob, use.mRNA, gi.type, f1.list=c("g1", "g2", "strata(types)"), f2.list=c("strata(sex)", "age", "strata(race)", "gii")) {

  # gi.type: one of "sl", "sdl", "dusr", "ddsr", "uusr", "udsr"

  if (use.mRNA){
    g1 <- prob$mRNA.norm[pair[1],]  
    g2 <- prob$mRNA.norm[pair[2],]
    f1 <- prob$mRNAq2[pair[1],]  
    f2 <- prob$mRNAq2[pair[2],]  
  } else {
    g1 <- prob$scna.norm[pair[1],]  
    g2 <- prob$scna.norm[pair[2],]  
    f1 <- prob$scnaq2[pair[1],]  
    f2 <- prob$scnaq2[pair[2],]  
  }

  dt1 <- cbind(prob$surv.strata, cbind(g1, g2))
  ll1 <- coxph.robust(f1.list=f1.list, f2.list=f2.list, data1=dt1)

  cov <- switch(gi.type,
                sdl=,
                dusr=ifelse(f2==0 & f1==2, 1, 0),
                sl=,
                ddsr=ifelse(f2==0 & f1==0, 1, 0),
                uusr=ifelse(f2==2 & f1==2, 1, 0),
                udsr=ifelse(f2==2 & f1==0, 1, 0))
  dt1$cov <- qnorm.array(cov)
  cox.out <- coxph.robust(f1.list=c(f1.list, "cov"), f2.list=f2.list, data1=dt1) 
  ll3 <- cox.out[6]
  bb  <- cox.out[1:5]
  uu <- c(bb, ll3-ll1)

  if (gi.type %in% c("dusr","ddsr","uusr","udsr")) {
    cov <- switch(gi.type,
                  dusr=ifelse(f2==0 & f1<=1, 1, 0),
                  ddsr=ifelse(f2==0 & f1>=1, 1, 0),
                  uusr=ifelse(f2==2 & f1<=1, 1, 0),
                  udsr=ifelse(f2==2 & f1>=1, 1, 0))
    dt1$cov <- qnorm.array(cov)
    cox.out <- coxph.robust(f1.list=c(f1.list, "cov"), f2.list=f2.list, data1=dt1) 
    ll2 <- cox.out[6]
    aa  <- cox.out[1:5]
    uu <- c(aa, ll2-ll1, uu)
  }
  
  return(uu)  

}

clin.screens <- function(queries, query.type="v", gi.type, use.mRNA, f1.list=c("g1", "g2", "strata(types)"), f2.list=c("strata(sex)", "age", "strata(race)", "gii")) {
  
  if (is.matrix(queries)) {
    sr <- queries
  } else if (is.vector(queries)) {
    queries.ind <- match(queries, genes)
    queries.ind <- queries.ind[!is.na(queries.ind)]
    tot.genes <- length(genes)
    n.queries <- length(queries.ind)
    if (query.type=="v") {
      sr <- cbind(ri=rep(1:tot.genes, each=n.queries), vi=rep(queries.ind, tot.genes))
    } else if (query.type=="r") {
      sr <- cbind(ri=rep(queries.ind, each=tot.genes), vi=rep(1:tot.genes, n.queries))
    }
  }

  res <- mclapply(1:nrow(sr), function(tt) {
    res <- clin.screen(sr[tt,], pancancer, use.mRNA=use.mRNA, gi.type="dusr", f1.list=f1.list, f2.list=f2.list)
    names(res) <- NULL
    gc()
    return(res)
  }, mc.cores=ncores)

  res <- as.data.table(do.call(rbind, res))
  res[, c("rn","vn","ri","vi"):=list(genes[sr[,1]], genes[sr[,2]], sr[,1], sr[,2])]
  setcolorder(res, c("rn","vn","ri","vi",names(res)[1:(length(res)-4)]))

  return(res)
}

cox.pair.du.strata.controlled = function(pair,prob, use.mRNA=F, f2.list = NULL, f1.list=NULL)
{
  if(use.mRNA){
    g1 = prob$mRNA.norm[pair[1],]  
    g2 = prob$mRNA.norm[pair[2],]
    f1 = prob$mRNAq2[pair[1],]  
    f2 = prob$mRNAq2[pair[2],]  
  }else{
    g1 = prob$scna.norm[pair[1],]  
    g2 = prob$scna.norm[pair[2],]  
    f1 = prob$scnaq2[pair[1],]  
    f2 = prob$scnaq2[pair[2],]  
  }
  surv.strata = prob$surv.strata
  if(is.null(f2.list)) f2.list = c(  "strata(sex)",  "age", "strata(race)") 
  if(is.null(f1.list)) f1.list =  "strata(types)"
  dt1 = cbind(surv.strata, cbind(g1 , g2))
  cntrl.out = coxph.robust(f1.list = c("g1",  "g2",  f1.list), f2.list = f2.list,data1=dt1)
  ll1 = cntrl.out
  cov = ifelse(f2 == 0 & f1<=1, 1, 0 )
  dt1$cov = qnorm.array(cov)
  cox.out = coxph.robust(f1.list = c("g1",  "g2",  f1.list, "cov"), f2.list = f2.list,data1=dt1) 
  ll2 = cox.out[6]
  aa  = cox.out[1:5]
  cov = ifelse(f2 == 0 & f1==2, 1, 0 )
  dt1$cov = qnorm.array(cov)
  cox.out = coxph.robust(f1.list = c("g1",  "g2",  f1.list, "cov"), f2.list = f2.list,data1=dt1) 
  ll3 = cox.out[6]
  bb  = cox.out[1:5]
  uu = c(aa, ll2 - ll1, bb, ll3-ll1)
  
  dt1 = surv.strata
  cov = ifelse(f2 == 0 & f1<=1, 1, 0 )
  dt1$cov = qnorm.array(cov)
  cox.out = coxph.robust(f1.list = c(  f1.list, "cov"), f2.list = f2.list,data1=dt1)  
  aa  = cox.out[1:5]
  cov = ifelse(f2 == 0 & f1==2, 1, 0 )
  dt1$cov = qnorm.array(cov)
  cox.out = coxph.robust(f1.list = c(  f1.list, "cov"), f2.list = f2.list,data1=dt1)  
  bb  = cox.out[1:5]
  c(uu, aa,  bb)
  
}

cox.pair.dd.strata.controlled = function(pair,prob, use.mRNA=F, f2.list = NULL, f1.list=NULL)
{
  if(use.mRNA){
    g1 = prob$mRNA.norm[pair[1],]  
    g2 = prob$mRNA.norm[pair[2],]
    f1 = prob$mRNAq2[pair[1],]  
    f2 = prob$mRNAq2[pair[2],]  
  }else{
    g1 = prob$scna.norm[pair[1],]  
    g2 = prob$scna.norm[pair[2],]  
    f1 = prob$scnaq2[pair[1],]  
    f2 = prob$scnaq2[pair[2],]  
  }
  surv.strata = prob$surv.strata
  if(is.null(f2.list)) f2.list = c(  "strata(sex)",  "age", "strata(race)") 
  if(is.null(f1.list)) f1.list =  "strata(types)"
  dt1 = cbind(surv.strata, cbind(g1 , g2))
  cntrl.out = coxph.robust(f1.list = c("g1",  "g2",  f1.list), f2.list = f2.list,data1=dt1)
  ll1 = cntrl.out
  cov = ifelse(f2 == 0 & f1>=1, 1, 0 )
  dt1$cov = qnorm.array(cov)
  cox.out = coxph.robust(f1.list = c("g1",  "g2",  f1.list, "cov"), f2.list = f2.list,data1=dt1) 
  ll2 = cox.out[6]
  aa  = cox.out[1:5]
  cov = ifelse(f2 == 0 & f1==0, 1, 0 )
  dt1$cov = qnorm.array(cov)
  cox.out = coxph.robust(f1.list = c("g1",  "g2",  f1.list, "cov"), f2.list = f2.list,data1=dt1) 
  ll3 = cox.out[6]
  bb  = cox.out[1:5]
  uu = c(aa, ll2 - ll1, bb, ll3-ll1)
  
  dt1 = surv.strata
  cov = ifelse(f2 == 0 & f1>=1, 1, 0 )
  dt1$cov = qnorm.array(cov)
  cox.out = coxph.robust(f1.list = c(  f1.list, "cov"), f2.list = f2.list,data1=dt1)  
  aa  = cox.out[1:5]
  cov = ifelse(f2 == 0 & f1==0, 1, 0 )
  dt1$cov = qnorm.array(cov)
  cox.out = coxph.robust(f1.list = c(  f1.list, "cov"), f2.list = f2.list,data1=dt1)  
  bb  = cox.out[1:5]
  c(uu, aa,  bb)
  
}

cox.pair.sl.strata.controlled = function(pair,prob, use.mRNA=F, f2.list = NULL, f1.list = NULL)
{

  if(use.mRNA){
    g1 = prob$mRNA.norm[pair[1],]  
    g2 = prob$mRNA.norm[pair[2],]
    f1 = prob$mRNAq2[pair[1],]  
    f2 = prob$mRNAq2[pair[2],]  
  }else{
    g1 = prob$scna.norm[pair[1],]  
    g2 = prob$scna.norm[pair[2],]  
    f1 = prob$scnaq2[pair[1],]  
    f2 = prob$scnaq2[pair[2],]  
  }
  surv.strata = prob$surv.strata
  if(is.null(f1.list)) f1.list =  "strata(type)"
  if(is.null(f2.list)) f2.list = c(  "strata(sex)",  "age", "strata(race)") 
  dt1 = cbind(surv.strata, cbind(g1 , g2))
  cntrl.out = coxph.robust(f1.list = c("g1",  "g2",  f1.list), f2.list = f2.list,data1=dt1)
  ll1 = cntrl.out
  cov = ifelse(f2 == 0 & f1 == 0, 1, 0 )
  dt1$cov = qnorm.array(cov)
  cox.out = coxph.robust(f1.list = c("g1",  "g2",  f1.list, "cov"), f2.list = f2.list,data1=dt1) 
  ll2 = cox.out[6]
  aa  = cox.out[1:5]
  
  uu = c(aa, ll2 - ll1)
  #dt1 = surv.strata
  #cov = ifelse(f2 == 0 & f1==0, 1, 0 )
  #dt1$cov = qnorm.array(cov)
  #cox.out = coxph.robust(f1.list = c(  f1.list, "cov"), f2.list = f2.list,data1=dt1)  
  #aa  = cox.out[1:5]
  #c(uu, aa)
  
} 


cox.pair.sdl.strata.controlled = function(pair,prob, use.mRNA=F, f2.list = NULL, f1.list = NULL)
{

  if(use.mRNA){
    g1 = prob$mRNA.norm[pair[1],]  
    g2 = prob$mRNA.norm[pair[2],]
    f1 = prob$mRNAq2[pair[1],]  
    f2 = prob$mRNAq2[pair[2],]  
  }else{
    g1 = prob$scna.norm[pair[1],]  
    g2 = prob$scna.norm[pair[2],]  
    f1 = prob$scnaq2[pair[1],]  
    f2 = prob$scnaq2[pair[2],]  
  }
  surv.strata = prob$surv.strata
  if(is.null(f1.list)) f1.list =  "strata(type)"
  if(is.null(f2.list)) f2.list = c(  "strata(sex)",  "age", "strata(race)") 
  dt1 = cbind(surv.strata, cbind(g1 , g2))
  cntrl.out = coxph.robust(f1.list = c("g1",  "g2",  f1.list), f2.list = f2.list,data1=dt1)
  ll1 = cntrl.out
  cov = ifelse(f2 == 0 & f1 == 2, 1, 0 )
  dt1$cov = qnorm.array(cov)
  cox.out = coxph.robust(f1.list = c("g1",  "g2",  f1.list, "cov"), f2.list = f2.list,data1=dt1) 
  ll2 = cox.out[6]
  aa  = cox.out[1:5]
  
  uu = c(aa, ll2 - ll1)
  #dt1 = surv.strata
  #cov = ifelse(f2 == 0 & f1==2, 1, 0 )
  #dt1$cov = qnorm.array(cov)
  #cox.out = coxph.robust(f1.list = c(  f1.list, "cov"), f2.list = f2.list,data1=dt1)  
  #aa  = cox.out[1:5]
  #c(uu, aa)
  
} 
