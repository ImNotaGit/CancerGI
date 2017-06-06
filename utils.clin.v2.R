library(data.table)
library(parallel)
library(survival)
ncores=23L

load("/cbcb/project2-scratch/kycheng/GI/tcga.genes.RData") # genes
load("/cbcb/project2-scratch/kycheng/GI/tcga.norm.RData") # scna.norm, mRNA.norm
load("/cbcb/project2-scratch/kycheng/GI/tcga.q2.RData") # scnaq2, mRNAq2
load("/cbcb/project2-scratch/kycheng/GI/tcga.surv.strata.RData") # pancancer.surv.strata

inv.norm <- function(vec) {
  # a function to normalize a vector to normal distribution
  vec.rank <- rank(vec, na.last="keep", ties.method="average")
  res <- vec.rank / (sum(!is.na(vec.rank)) + 1)
  qnorm(res)
}

make.cox.formula <- function(vars, data="pancancer.surv.strata") {
  library(stringr)
  
  tmpl <- grepl("strata", vars)
  tmp <- vars[tmpl]
  tmp <- str_replace(tmp, "strata\\(", paste0("strata\\(",data,"$"))
  vars[tmpl] <- tmp
  
  tmpl <- !vsl & !(vars %in% c("cov","r","v"))
  tmp <- vars[tmpl]
  tmp <- paste0(data, "$", tmp)
  vars[tmpl] <- tmp

  as.formula(paste(sprintf("Surv(%s$time, %s$status) ~",data,data), paste(vars, collapse="+")))
}

clin.screen <- function(queries, query.type="v", gi.type, use.mRNA, vars=c("r", "v", "strata(types)"), vars.xtra=c("age", "strata(sex)", "strata(race)", "gii")) {
  # queries should be a vector of gene symbols, in this case query.type should be "v" so that the queries are vulnerable genes or "r" so that the queries are rescuers
  # or, queries can be a data.table of two columns containing pairs of genes, the column names should be vi and ri (gene indeces), or vn and rn (gene symbols)
  # gi.type: one of "sl", "sdl", "dusr", "ddsr", "uusr", "udsr"

  if (is.data.table(queries)) {
    # if queries is a data.table, assume it contains gene pairs
    if ("ri" %in% names(queries) & "vi" %in% names(queries)) {
      # assume gene indeces are provided in this case, and assume they are in the range of genes
      pairs <- transpose(queries[, .(ri, vi)])
      cols.add <- cbind(queries[, .(rn=genes[ri], vn=genes[vi])], queries)
    } else if ("rn" %in% names(queries) & "vn" %in% names(queries)) {
      # assume gene symbols are provided in this case
      cols.add <- copy(queries)
      cols.add[, c("ri", "vi"):=list(match(rn,genes), match(vn,genes))]
      cols.add <- cols.add[!is.na(ri) & !is.na(vi)]
      pairs <- transpose(cols.add[, .(ri, vi)])
    }
  } else if (is.vector(queries)) {
    # else if queries is a vector, assume it's a vector of gene symbols
    qinds <- match(queries, genes)
    qinds <- qinds[!is.na(qinds)]
    nq <- length(qinds)
    ntot <- length(genes)
    if (query.type=="v") {
      ri <- rep(1:ntot, each=nq)
      vi <- rep(qinds, ntot)
    } else if (query.type=="r") {
      ri <- rep(qinds, each=ntot)
      vi <- rep(1:ntot, nq)
    } else stop("query.type not given as appropriate.\n")
    pairs <- as.data.table(rbind(ri, vi))
    cols.add <- data.table(rn=genes[ri], vn=genes[vi], ri=ri, vi=vi)
  }

  formula.try1 <- make.cox.formula(c(vars, vars.xtra))
  formula.cov.try1 <- make.cox.formula(c("cov", vars, vars.xtra))
  formula.try2 <- make.cox.formula(vars)
  formula.cov.try2 <- make.cox.formula(c("cov", vars))

  cox <- function(pair) {
      r.ind <- pair[1]
      v.ind <- pair[2]
    if (use.mRNA){
      r <- mRNA.norm[r.ind,]
      v <- mRNA.norm[v.ind,]
      rq <- mRNAq2[r.ind,]
      vq <- mRNAq2[v.ind,]
    } else {
      r <- scna.norm[r.ind,]
      v <- scna.norm[v.ind,]
      rq <- scnaq2[r.ind,]
      vq <- scnaq2[v.ind,]
    }
  
    lik0 <- tryCatch(coxph(formula.try1), error=function(e) coxph(formula.try2))$loglik[2]
  
    cov <- switch(gi.type,
                  sdl=,
                  dusr=ifelse(vq==0 & rq==2, 1, 0),
                  sl=,
                  ddsr=ifelse(vq==0 & rq==0, 1, 0),
                  uusr=ifelse(vq==2 & rq==2, 1, 0),
                  udsr=ifelse(vq==2 & rq==0, 1, 0))
    cov <- inv.norm(cov)
    tmp <- tryCatch(coxph(formula.cov.try1), error=function(e) coxph(formula.cov.try2))
    tmp <- c(summary(tmp)$coefficients["cov",], tmp$loglik[2])
    res <- c(tmp[1:5], 1-pchisq(2*(lik0-tmp[6]),1))
  
    if (gi.type %in% c("dusr","ddsr","uusr","udsr")) {
      cov <- switch(gi.type,
                    dusr=ifelse(vq==0 & rq<=1, 1, 0),
                    ddsr=ifelse(vq==0 & rq>=1, 1, 0),
                    uusr=ifelse(vq==2 & rq<=1, 1, 0),
                    udsr=ifelse(vq==2 & rq>=1, 1, 0))
      cov <- inv.norm(cov)
      tmp <- tryCatch(coxph(formula.cov.try1), error=function(e) coxph(formula.cov.try2))
      tmp <- c(summary(tmp)$coefficients["cov",], tmp$loglik[2])
      res <- c(tmp[1:5], 1-pchisq(2*(lik0-tmp[6]),1), res)
    }
    
    return(res)
  }

  p <- as.data.table(do.call(rbind, mclapply(pairs, cox, mc.cores=ncores)))
  cbind(cols.add, p)

}
