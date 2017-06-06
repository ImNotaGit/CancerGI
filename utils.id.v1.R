get.rawp.gsof <- function(queries, query.type="v", gi.type="dusr", data.dir="/cbcb/project2-scratch/kycheng/GI/gsof") {
  
  # get pre-calculated gSOF raw p values for query-by-all gene pairs.
  # queries: a vector of query genes (as gene symbols)
  # query.type: query genes as vulnerable partners ("v") or rescuers ("r"); for sl, these make no difference except for the resulting layout of data; for sdl, use query.type in the same way as DUSR: "v" means the query genes are the "down genes", "r" means the query genes are the "up genes".
  # gi.type: type of gene interaction to identify, e.g. "sl", "sdl", "dusr", etc.
  # data.dir: the dir containing pre-calculated gSOFs and gene data
  
  library(data.table)
  library(stringr)
  
  # remove trailing "/" in data.dir, and if empty, set to "."
  data.dir <- str_replace(data.dir, "/$", "")
  if (data.dir=="") data.dir <- "."
  
  # all genes
  load(file.path(data.dir, "genes.RData"))
  tot.genes <- length(genes)
  
  # query genes
  queries.ind <- match(queries, genes)
  queries.ind <- queries.ind[!is.na(queries.ind)]
  queries.name <- genes[queries.ind]
  n.queries <- length(queries.ind)
  if (n.queries==0) stop("No match for any of the query genes.\n")
  # indeces in the pre-calculated gsof data for query-by-all gene pairs
  if (query.type=="v") {
    queries.pairs.ind <- as.vector(sapply(0:(tot.genes-1), function(x) x*tot.genes+queries.ind))
    v.inds <- rep(queries.ind, tot.genes)
    r.inds <- rep(1:tot.genes, each=n.queries)
  } else if (query.type=="r") {
    queries.pairs.ind <- as.vector(sapply(queries.ind-1, function(x) x*tot.genes+(1:tot.genes)))
    v.inds <- rep(1:tot.genes, n.queries)
    r.inds <- rep(queries.ind, each=tot.genes)
  }
  v.names <- genes[v.inds]
  r.names <- genes[r.inds]
  
  # load pre-calculated gsof data
  file.inds <- switch(gi.type,
                      sl=0,
                      sdl=2,
                      dusr=0:2,
                      ddsr=0:2,
                      uusr=6:8,
                      udsr=6:8)
  data.files <- file.path(data.dir, outer(c("scna","mrna"), paste0(".gsof", file.inds, ".RData"), paste0))
  gsof.summ <- as.data.table(lapply(data.files, function(fn) {
    obj.name <- load(fn)
    obj <- get(obj.name)
    obj[queries.pairs.ind]
  }))
  
  if (gi.type %in% c("sl","sdl")) {
    # set column names
    setnames(gsof.summ, c("scna","mrna"))
  } else if (gi.type %in% c("dusr","ddsr","uusr","udsr")) {
    # set column names
    setnames(gsof.summ, outer(c("scna.r.","mrna.r."), c("lo","med","hi"), paste0))
    # reverse the direction of specific columns of p values based on the wanted type of GI
    if (gi.type %in% c("dusr","uusr")) {
      gsof.summ[, c("scna.r.hi", "mrna.r.hi"):=list(1-scna.r.hi, 1-mrna.r.hi)]
    } else if (gi.type %in% c("ddsr","udsr")) {
      gsof.summ[, c("scna.r.lo", "mrna.r.lo"):=list(1-scna.r.lo, 1-mrna.r.lo)]
    }
  } else stop("gi.type not given as appropriate.\n")

  # add the indeces and symbols for the gene pairs
  gsof.summ[, c("rn","vn","ri","vi"):=list(r.names, v.names, r.inds, v.inds)]
  setcolorder(gsof.summ, c("rn","vn","ri","vi",names(gsof.summ)[1:(length(gsof.summ)-4)]))

  # keep recorde of the gi type as an attribute. This can be easily lost, though.
  attr(gsof.summ, "gi.type") <- gi.type
  return(gsof.summ)
  
}


id.gi.gsof <- function(rawp, queries=NULL, query.type=NULL, gi.type="sr", fdr.thr=0.05) {
  
  # Identify gene interactions based on pre-calculated gSOFs.
  # rawp: a data.table containing raw p values returned from get.rawp.gsof()
  # queries: a vector of query genes (as gene symbols); if NULL, all records in rawp are used
  # query.type: query genes as vulnerable partners ("v") or rescuers ("r"); for sdl, this is similar to DUSR; for sl, this should correspond to (i.e, be the same as) the choice in get.rawp.gsof().
  # gi.type: "sr" or "sl"
  # fdr.thr: FDR threshold
  
  library(data.table)
  
  if (gi.type=="sl") {
    tmp <- copy(rawp)
  } else if (gi.type=="sr") {
    tmp <- rawp[, .(rn, vn, ri, vi, scna=pmin(scna.r.lo,scna.r.hi), mrna=pmin(mrna.r.lo,mrna.r.hi))]
  } else stop("gi.type not given as appropriate.\n")

  if (!is.null(queries)) {
    if (query.type=="v") {
      tmp <- tmp[vn %in% queries]
    } else if (query.type=="r") {
      tmp <- tmp[rn %in% queries]
    } else stop("query.type not given as appropriate.\n")
  }
  
  tmp[, prod:=scna*mrna]
  tmp[, score:=ifelse(scna<0.05 & mrna<0.05, prod, 1)]
  # FDR correction
  tmp[, score:=p.adjust(score, method="BH")]
  tmp <- tmp[score<fdr.thr, .(rn, vn, ri, vi, score)]
  # return result
  if (is.null(queries)) {
    return(tmp[order(score)])
  } else if (query.type=="v") {
    return(tmp[order(vi, score)])
  } else if (query.type=="r") {
    return(tmp[order(ri, score)])
  }
  
}


id.gi.clin <- function(rawp.scna, rawp.mrna, queries=NULL, query.type=NULL, gi.type="sr", fdr.thr=0.05) {
  
  # Identify gene interactions based on pre-calculated p values for clinical screen
  # rawp.scna, rawp.mrna: data.tables containing raw p values for scna and mRNA returned from the standard pipeline
  # queries: a vector of query genes (as gene symbols); if NULL, all records in rawp are used
  # query.type: query genes as vulnerable partners ("v") or rescuers ("r"); for sdl, this is similar to DUSR; for sl, this should be consistent with the choice in the previous screening steps
  # gi.type: "sl" or "sr"
  # fdr.thr: FDR threshold
  
  library(data.table)

  if (is.null(queries)) {
    tmp.scna <- copy(rawp.scna)
    tmp.mrna <- copy(rawp.mrna)
  } else {
    if (query.type=="v") {
      tmp.scna <- rawp.scna[vn %in% queries]
      tmp.mrna <- rawp.mrna[vn %in% queries]
    } else if (query.type=="r") {
      tmp.scna <- rawp.scna[rn %in% queries]
      tmp.mrna <- rawp.mrna[rn %in% queries]
    } else stop("query.type not given as appropriate.\n")
  }

  if (gi.type=="sr") {
    # V1, V7: coefficients (beta) of non-rescued bin (e.g. D-not.U in DUSR) and rescued bin (e.g. D-U in DUSR); V5, V11: p values of the two bins
    tmp.scna <- tmp.scna[, .(rn,vn,ri,vi,V1,V7,V5,V11)]
    tmp.mrna <- tmp.mrna[, .(rn,vn,ri,vi,V1,V7,V5,V11)]
    tmp.scna[, score.scna:=ifelse(V1<0 & V7>0, V5*V11, 1)]
    tmp.mrna[, score.mrna:=ifelse(V1<0 & V7>0, V5*V11, 1)]
  } else if (gi.type=="sl") {
    # V1: coefficients (beta) for the sl bin; V5: p values for the sl bin
    tmp.scna <- tmp.scna[, .(rn,vn,ri,vi,V1,V5)]
    tmp.mrna <- tmp.mrna[, .(rn,vn,ri,vi,V1,V5)] 
    tmp.scna[, score.scna:=ifelse(V1<0, V5, 1)]
    tmp.mrna[, score.mrna:=ifelse(V1<0, V5, 1)]
  } else stop("gi.type not given as appropriate.\n")

  # fdr correction
  tmp.scna[, score.scna:=p.adjust(score.scna, method="BH")]
  tmp.scna <- tmp.scna[score.scna<fdr.thr, .(rn, vn, ri, vi, score.scna)]
  tmp.mrna[, score.mrna:=p.adjust(score.mrna, method="BH")]
  tmp.mrna <- tmp.mrna[score.mrna<fdr.thr, .(ri, vi, score.mrna)]

  # result
  res <- tmp.scna[tmp.mrna, .(rn, vn, ri, vi, score.scna, score.mrna), on=c("ri","vi"), nomatch=0]
  # return
  if (is.null(queries)) {
    return(res[order(pmin(score.scna, score.mrna))])
  } else if (query.type=="v") {
    return(res[order(vi, pmin(score.scna, score.mrna))])
  } else if (query.type=="r") {
    return(res[order(ri, pmin(score.scna, score.mrna))])
  }

}

id.gi.clin.nofdr <- function(rawp.scna, rawp.mrna, queries=NULL, query.type=NULL, fdr.thr=0.05) {
  
  # Identify gene interactions (for now, SR only) based on pre-calculated p values for clinical screen
  # rawp.scna, rawp.mrna: data.tables containing raw p values for scna and mRNA returned from the standard pipeline
  # queries: a vector of query genes (as gene symbols); if NULL, all records in rawp are used
  # query.type: query genes as vulnerable partners ("v") or rescuers ("r")
  # fdr.thr: FDR threshold
  
  library(data.table)

  # V1, V7: coefficients (beta) of non-rescued bin (e.g. D-not.U in DUSR) and rescued bin (e.g. D-U in DUSR); V5, V11: p values of the two bins
  tmp.scna <- rawp.scna[, .(rn,vn,ri,vi,V1,V7,V5,V11)]
  tmp.mrna <- rawp.mrna[, .(rn,vn,ri,vi,V1,V7,V5,V11)]
  if (!is.null(queries)) {
    if (query.type=="v") {
      tmp.scna <- tmp.scna[vn %in% queries]
      tmp.mrna <- tmp.mrna[vn %in% queries]
    } else if (query.type=="r") {
      tmp.scna <- tmp.scna[rn %in% queries]
      tmp.mrna <- tmp.mrna[rn %in% queries]
    } else stop("query.type not given as appropriate.\n")
  }
  # fdr correction
  tmp.scna[, score.scna:=ifelse(V1<0 & V7>0, V5*V11, 1)]
#  tmp.scna[, score.scna:=p.adjust(score.scna, method="BH")]
  tmp.scna <- tmp.scna[score.scna<fdr.thr, .(rn, vn, ri, vi, score.scna)]
  tmp.mrna[, score.mrna:=ifelse(V1<0 & V7>0, V5*V11, 1)]
#  tmp.mrna[, score.mrna:=p.adjust(score.mrna, method="BH")]
  tmp.mrna <- tmp.mrna[score.mrna<fdr.thr, .(ri, vi, score.mrna)]
  # result
  res <- tmp.scna[tmp.mrna, .(rn, vn, ri, vi, score.scna, score.mrna), on=c("ri","vi"), nomatch=0]
  # return
  if (is.null(queries)) {
    return(res[order(pmin(score.scna, score.mrna))])
  } else if (query.type=="v") {
    return(res[order(vi, pmin(score.scna, score.mrna))])
  } else if (query.type=="r") {
    return(res[order(ri, pmin(score.scna, score.mrna))])
  }

}

id.gi.ess.du <- function(rawp, queries=NULL, query.type=NULL, fdr.thr=0.05) {

  # Identify gene interactions (this version of the function is for DUSR) based on pre-calculated p values from gene essentiality screens
  # rawp: a list of data.tables containing raw p values from multiple screening datasets
  # queries: a vector of query genes (as gene symbols); if NULL, all records in rawp are used
  # query.type: query genes as vulnerable partners ("v") or rescuers ("r")
  # fdr.thr: FDR threshold
  
  library(data.table)

  if (is.null(queries)) {
    tmp <- rbindlist(rawp, idcol="screen.id")[, .(screen.id, rn,vn,ri,vi, pmrna.v=V2, pscna.v=V4, pmrna.r=V1, pscna.r=V3)]
  } else {
    if (query.type=="v") {
      tmp <- rbindlist(rawp, idcol="screen.id")[vn %in% queries, .(screen.id, rn,vn,ri,vi, pmrna.v=V2, pscna.v=V4, pmrna.r=V1, pscna.r=V3)]
    } else if (query.type=="r") {
      tmp <- rbindlist(rawp, idcol="screen.id")[rn %in% queries, .(screen.id, rn,vn,ri,vi, pmrna.v=V2, pscna.v=V4, pmrna.r=V1, pscna.r=V3)]
    } else stop("query.type not given as appropriate.\n")
  }
  # fdr correction
  tmp[, c("padj.mrna.v","padj.scna.v","padj.mrna.r","padj.scna.r"):=lapply(.SD, p.adjust, method="BH"), .SDcols=c("pmrna.v","pscna.v","pmrna.r","pscna.r"), by=screen.id]
  # positive SRs: those with all four padj's < threshold in any one or more screening datasets
  res <- unique(tmp[pmin(padj.mrna.v, padj.scna.v, padj.mrna.r, padj.scna.r)<fdr.thr, .(rn,vn,ri,vi)])

  #res <- dcast(res, rn+vn+ri+vi~screen.id, sep=".", value.var=c("padj.mrna.v","padj.scna.v","padj.mrna.r","padj.scna.r")) # this will include the various padj's in all the screens. for this to work, screen.id need to be included in res in the line above.

  # return
  if (is.null(queries)) {
    return(res)
  } else if (query.type=="v") {
    return(res[order(vi)])
  } else if (query.type=="r") {
    return(res[order(ri)])
  }
  
}

id.gi.ess.sl <- id.gi.ess.du

id.gi.ess.dd <- function(rawp, queries=NULL, query.type=NULL, fdr.thr=0.05) {
  
  # Identify gene interactions (this version of the function is for DDSR) based on pre-calculated p values from gene essentiality screens
  # rawp: a list of data.tables containing raw p values from multiple screening datasets
  # queries: a vector of query genes (as gene symbols); if NULL, all records in rawp are used
  # query.type: query genes as vulnerable partners ("v") or rescuers ("r")
  # fdr.thr: FDR threshold
  
  library(data.table)

  if (is.null(queries)) {
    tmp <- rbindlist(rawp, idcol="screen.id")[, .(screen.id, rn,vn,ri,vi, pmrna.r=V5, pscna.r=V7)]
  } else {
    if (query.type=="v") {
      tmp <- rbindlist(rawp, idcol="screen.id")[vn %in% queries, .(screen.id, rn,vn,ri,vi, pmrna.r=V5, pscna.r=V7)]
    } else if (query.type=="r") {
      tmp <- rbindlist(rawp, idcol="screen.id")[rn %in% queries, .(screen.id, rn,vn,ri,vi, pmrna.r=V5, pscna.r=V7)]
    } else stop("query.type not given as appropriate.\n")
  }
  # fdr correction
  tmp[, c("padj.mrna.r","padj.scna.r"):=lapply(.SD, p.adjust, method="BH"), .SDcols=c("pmrna.r","pscna.r"), by=screen.id]
  # positive SRs: those with all four padj's < threshold in any one or more screening datasets
  res <- unique(tmp[pmin(padj.mrna.r, padj.scna.r)<fdr.thr, .(rn,vn,ri,vi)])

  #res <- dcast(res, rn+vn+ri+vi~screen.id, sep=".", value.var=c("padj.mrna.v","padj.scna.v","padj.mrna.r","padj.scna.r")) # this will include the various padj's in all the screens. for this to work, screen.id need to be included in res in the line above.

  # return
  if (is.null(queries)) {
    return(res)
  } else if (query.type=="v") {
    return(res[order(vi)])
  } else if (query.type=="r") {
    return(res[order(ri)])
  }
  
}

id.gi.ess.uu <- function(rawp, queries=NULL, query.type=NULL, gi.type="dusr", fdr.thr=0.05) {

  # Identify gene interactions (this version of the function is for UUSR) based on pre-calculated p values from gene essentiality screens
  # rawp: a list of data.tables containing raw p values from multiple screening datasets
  # queries: a vector of query genes (as gene symbols); if NULL, all records in rawp are used
  # query.type: query genes as vulnerable partners ("v") or rescuers ("r")
  # fdr.thr: FDR threshold
  
  library(data.table)

  if (is.null(queries)) {
    tmp <- rbindlist(rawp, idcol="screen.id")[, .(screen.id, rn,vn,ri,vi, pmrna.v=V6, pscna.v=V8)]
  } else {
    if (query.type=="v") {
      tmp <- rbindlist(rawp, idcol="screen.id")[vn %in% queries, .(screen.id, rn,vn,ri,vi, pmrna.v=V6, pscna.v=V8)]
    } else if (query.type=="r") {
      tmp <- rbindlist(rawp, idcol="screen.id")[rn %in% queries, .(screen.id, rn,vn,ri,vi, pmrna.v=V6, pscna.v=V8)]
    } else stop("query.type not given as appropriate.\n")
  }
  # fdr correction
  tmp[, c("padj.mrna.v","padj.scna.v"):=lapply(.SD, p.adjust, method="BH"), .SDcols=c("pmrna.v","pscna.v"), by=screen.id]
  # positive SRs: those with all four padj's < threshold in any one or more screening datasets
  res <- unique(tmp[pmin(padj.mrna.v, padj.scna.v)<fdr.thr, .(rn,vn,ri,vi)])

  #res <- dcast(res, rn+vn+ri+vi~screen.id, sep=".", value.var=c("padj.mrna.v","padj.scna.v","padj.mrna.r","padj.scna.r")) # this will include the various padj's in all the screens. for this to work, screen.id need to be included in res in the line above.

  # return
  if (is.null(queries)) {
    return(res)
  } else if (query.type=="v") {
    return(res[order(vi)])
  } else if (query.type=="r") {
    return(res[order(ri)])
  }
  
}