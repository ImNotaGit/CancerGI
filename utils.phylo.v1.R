load("/cbcb/project2-scratch/kycheng/GI/tcga.genes.RData") # genes
# The phylogenetic profile is downloaded from Yuval Tabach et al. Mol Syst Biol. (2013), Supplementary Table 1
load("/cbcb/project2-scratch/kycheng/GI/yuval.phylogenetic.profile.RData") # phylo
# the feature weights are determined based on the phylogenetic tree (Ensembl database: http://useast.ensembl.org/index.html)
load("/cbcb/project2-scratch/kycheng/GI/feature.weight.RData") # feature.weight

# calculate the phylogenetic distance of random gene pairs for false discovery correction, run only once and save the result
#sr.rand <- data.table(rescuer=sample(genes, 100000, replace=TRUE), vulnerable=sample(genes, 100000, replace=TRUE))
#sr.rand <- sr.rand[rescuer!=vulnerable]
#pp.rand <- phylo.profile(sr.rand)
#save(pp.rand, file="phylo.rand.RData")
load("/cbcb/project2-scratch/kycheng/GI/phylo.rand.RData") # pp.rand

phylo.profile = function(sr.gene.all){
  load("/cbcb/project2-scratch/jooslee/srescues/pancancer/cox.du/yuval.phylogenetic.profile.RData")
  load("/cbcb/project2-scratch/jooslee/srescues/pancancer/cox.du/feature.weight.RData")
  sr.gene1 = sr.gene.all
  sr.phylo =  cbind(match(sr.gene1$rescuer, phylo$genes), match(sr.gene1$vulnerable, phylo$genes))
  featureMat = (phylo[sr.phylo[,1],-(1:3)] - phylo[sr.phylo[,2],-(1:3)])^2
  featureMat %*% t(feature.weight)
}