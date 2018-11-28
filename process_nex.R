library(lattice)
library(ape)
library(geiger)
library(phangorn)
library(phytools)
library(phyloclim)

#turk <- read.nexus.data(file = "data/Turkic12langs.nex")
turk <- read.nexus.data(file = "data/denominations.nex")
turk.phydata <- phyDat(as.data.frame(turk), format="NEXUS", type="USER", levels=c("0", "1"))
#print(as.data.frame(turk))
dm1 <- dist.ml(turk.phydata)
dm2 <- dist(t(as.matrix(as.data.frame(lapply(turk, as.numeric)))), method="euclidean")

