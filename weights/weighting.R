source("weights/weighting_scheme.R")

# Input file: phylogenetic tree with edge.lenghts
phylotree <- read.newick(file = "weights/ensembl_validated.nwk")
tips <- timetree[["tip.label"]]
timetree <- as_tibble(phylotree)

# df <- timetree[order(timetree$branch.length, decreasing = T),]

################## Prunning the species tree ############################

## Removing individual species
# tips <- tips[tips != "Pan_paniscus"]
# tips <- tips[tips != "Macaca_fascicularis"]
tips <- tips[tips != "Ornithorhynchus_anatinus"]
tips <- tips[tips != "Xenopus_tropicalis"]

tips <- c(tips, "Xenopus_tropicalis")

## Removing all but OHNOLOGSv2 vertebrates
# ssp = readLines("weights/vertebrates_2r.txt")
# tips <- tips[(tips %in% ssp)]

rem <- c("Papio_anubis", "Theropithecus_gelada", "Macaca_mulatta", 'Macaca_fascicularis', "Chlorocebus_sabaeus")
tips <- tips[!(tips %in% rem)]


#########################################################################

pesos <- calculate_weights(timetree = timetree, tips = tips)
sum(pesos)

library(ggtree)

ggtree(phylotree) + geom_tiplab(as_ylab=TRUE, color='firebrick')

