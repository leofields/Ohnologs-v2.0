
source("weights/weighting_scheme.R")
library(ggtree)
library(dplyr)
library(ggstance)



# Input file:
phylotree <- read.newick(file = "weights/ensembl_validated.nwk")
tips <- phylotree[["tip.label"]]
species <- tips
timetree <- as_tibble(phylotree)

# Add OHNOLOGS_V2
v2 = readLines("weights/vertebrates_2r.txt")
timetree$ohnologs_v2 <- !is.na(match(timetree$label, v2))
# sum(timetree$ohnologs_v2)

# Species NOT MATCHED:
# v2[is.na(match(v2, dt$label))]
# 
# [1] "Callithrix_jacchus"  "Gasterosteus_aculeatus"


# Add common names
load("weights/AllSpecies.RData")
ensembl_species <- ensembl_species[!duplicated(ensembl_species$scientific_name),]
ensembl_species$label <- ensembl_species$scientific_name
t <- left_join(timetree, ensembl_species, by = "label")
t <- t[, c("label", "common_name")]

# Calculating the weights
spw <- calculate_weights(timetree = timetree, tips = tips)
sum(spw)

# Ploting the phylo tree

dt <- as_tibble(spw, rownames = "label")
colnames(dt)[2] <- "weight"
dt$weight <- round(dt$weight, 4)
dt$ohonologs_v2 <- !is.na(match(dt$label, v2))
sum(dt$ohonologs_v2)

tt <- full_join(timetree, dt, by = "label")
tt <- left_join(tt, t, by = "label")
tt$species <- tt$label
tt$name <- ifelse(is.na(tt$common_name), tt$species, tt$common_name)
# tt$label <- paste0("[", tt$weight, "] ", tt$name)
tt$label <- tt$name
tt$ohnologs_v2 <- !is.na(match(tt$species, v2))

p <- ggtree(phylotree) + geom_tiplab(align = T)
# ggtree(as.phylo(tt)) + xlim(-.1, 700)
p + geom_facet(panel = "Weight", data = dt, geom = ggstance::geom_barh, 
                aes(x = weight, color = ohonologs_v2), 
                stat = "identity", width = .6) +
    theme_tree2(legend.position=c(.05, .75))

# tr <- as.treedata(tt)
# 
# p <- ggtree(tr) + geom_tiplab() + xlim(-.1, 1000)
# p + geom_facet(panel = "Weights", data = dt, geom = geom_text, mapping = aes(label = weight))



## Removing individual species and re-calculating
tips <- tips[tips != "Pan_paniscus"]
spw <- calculate_weights(timetree = timetree, tips = tips)
sum(spw)

tips <- tips[tips != "Pan_troglodytes"]
spw <- calculate_weights(timetree = timetree, tips = tips)
sum(spw)

tips <- tips[tips != "Macaca_fascicularis"]
spw <- calculate_weights(timetree = timetree, tips = tips)
sum(spw)


rem <- c("Papio_anubis", "Theropithecus_gelada", "Macaca_mulatta", 'Macaca_fascicularis', "Chlorocebus_sabaeus")
tips <- tips[!(tips %in% rem)]
spw <- calculate_weights(timetree = timetree, tips = tips)
sum(spw)

dt2 <- as_tibble(spw, rownames = "label")
colnames(dt2)[2] <- "weight"
dt2$weight <- round(dt2$weight, 4)

# Combined plot
p + geom_facet(panel = "Weights ALL", data = dt, geom = ggstance::geom_barh, 
               aes(x = weight), 
               stat = "identity", width = .6) +
    geom_facet(panel = "_PRIMATES", data = dt2, geom = ggstance::geom_barh, 
               aes(x = weight), 
               stat = "identity", width = .6) +
    theme_tree2(legend.position=c(.05, .75))





tips <- tips[tips != "Ornithorhynchus_anatinus"]
spw <- calculate_weights(timetree = timetree, tips = tips)
sum(spw)

tips <- tips[tips != "Xenopus_tropicalis"]
spw <- calculate_weights(timetree = timetree, tips = tips)
sum(spw)


tips <- tips[tips != "Erpetoichthys_calabaricus"]
spw <- calculate_weights(timetree = timetree, tips = tips)
sum(spw)

tips <- tips[tips != "Monodelphis_domestica"]
spw <- calculate_weights(timetree = timetree, tips = tips)
sum(spw)



###########################
## Sorting from greatest to smallest branch lenghts
###########################

sorted <- timetree %>% 
    arrange(desc(branch.length)) %>% 
    filter(nchar(label) > 3)

###########################
## Annotate internal nodes
###########################
total_sum <- sum(spw)
nw <- data.frame(simple=rep(NA, nrow(tt)), actual=rep(NA, nrow(tt)))
for (n in 1:nrow(tt)) {
    # root of the tree
    if (tt[n, c('node')] == tt[n, c('parent')]){
        nw[n,] <- c(total_sum, total_sum)
        next
    }
    # tips of the tree 
    if (!is.numeric(tt[n, c('label')])){
        x <- tt[n, c('label')]
        y <- species[species != x]
        z <- calculate_weights(timetree, y)
        nw[n, ] <- c(tt[n, c('weight')], total_sum - sum(z))
    } else {
        # internal nodes
        x <- offspring(tt, n) %>% 
            filter(is.na(as.numeric(label))) %>% 
            select(species, weight)
        y <- species[!(species %in% x$species)]
        z <- calculate_weights(timetree, y)
        nw[n, ] <- c(sum(x$weight), total_sum - sum(z))
    }
}
annotated_tree <- bind_cols(tt,nw)


# simple_sum <- function(t,n) {
#     
#     aux <- offspring(t, n)
#     if(nrow(aux) == 0){
#         return(t[t$label == n,]$weight)
#     }
#     res <- aux %>% 
#         filter(is.na(as.numeric(label))) %>% 
#         select(weight)
#     return(sum(res))
# }
# mutate(tt, nw_simple = simple_sum(tt,label))



tt[74, c('node')] == tt[74, c('parent')]

teste <- offspring(tt, 1)

teste <- data.frame(N=numeric(), alpha=numeric())

teste <- offspring(tt, 144) %>% 
    filter(is.na(as.numeric(label))) %>% 
    select(species, weight)
teste

sum(teste)

tt$node == tt$parent

is.na(as.numeric(tt$label))

