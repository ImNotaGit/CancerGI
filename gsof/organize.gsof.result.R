#!/opt/local/stow/R-3.3.1/bin/Rscript

setwd("/cbcb/project2-scratch/kycheng/cancer_SR_drug/gsof")

# splitting the 9 cells into 0:8
for (i in 1:4) {
    load(paste0("scna",i,".RData"))
for (j in 1:9) {
      x <- scna[,j]
    save(x, file=paste0("scna",i,j-1,".RData"))
        rm(x)
        gc()
          
}
  rm(scna)
    gc()
  load(paste0("mrna",i,".RData"))
  for (j in 1:9) {
        x <- mrna[,j]
      save(x, file=paste0("mrna",i,j-1,".RData"))
          rm(x)
          gc()
            
  }
    rm(mrna)
      gc()

}

# joining the 4 parts for each cell
for (i in 0:8) {
  p <- unlist(lapply(1:4, function(j) {
                           load(paste0("scna",j,i,".RData"))
                               x
                             
}))
    save(p, file=paste0("scna.gsof",i,".RData"))
    rm(p)
      gc()
    p <- unlist(lapply(1:4, function(j) {
                              load(paste0("mrna",j,i,".RData"))
                                  x
                                
}))
      save(p, file=paste0("mrna.gsof",i,".RData"))
      rm(p)
        gc()

}

