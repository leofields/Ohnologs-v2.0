# Weighting scheme for phylogenetically related sequences
#
# This function requires as input a phylogenetic tree with branch lenghts.
# One such tree can be built from TimeTree website (http://timetree.org/)
# and exported/saved as a Newick file.


library(tidytree)
library(treeio)

# # Input file: phylogenetic tree with edge.lenghts
# timetree <- read.newick(file = "weights/ensembl_validated.nwk")
# tips <- timetree[["tip.label"]]
# timetree <- as_tibble(timetree)

calculate_weights <- function(timetree, tips) {
    
    #################### Auxiliar functions #################################
    #########################################################################
    # returns the branch length given a node label
    get_length <- function(tree, label) {
        round(tree[tree$label == label,]$branch.length, 1)
    }
    
    # returns the time of divergence between two species in the tree
    get_divergence_time <-  function(x, t){
        a <- MRCA(t, x[1], x[2])
        b <- parent(t, x[1])
        res <- get_length(t, x[1])
        while(a$label != b$label){
            res <- res + get_length(t, b$label)
            b <- parent(t, b$label)
        }
        res
    }
    #########################################################################
    
    # Pairing the species BY COMBINATIONS without repetitions
    ls <- as.list(as.data.frame(combn(tips, 2)))
    
    names(ls) <- sapply(ls, function(x) {
        paste0(sort(x)[1], sort(x)[2])
    })
    
    # obtaining the divergence time for each pair in the list
    div_list <- lapply(ls, get_divergence_time, timetree)
    
    ## Building the final matrix BY PERMUTATIONS with repetitions
    # Expanding a grid with all pairs of species
    aux <- expand.grid(rep(list(tips), 2), stringsAsFactors = F)
    
    # Fill the grid with divergence times from the list or zero in case of same sp.
    div_times <- apply(aux, 1, function(x){
        if(x[1] != x[2]){
            idx <- paste0(sort(x)[1], sort(x)[2])
            return (div_list[[idx]])
        }
        return (0)
    })
    
    mx <- matrix(div_times, nrow = length(tips), dimnames = list(tips, tips))
    
    ## Calculation of the weights
    det(mx)
    m = 535 - mx 
    m = m * m
    m = m / (535 * 535)
    
    inv = solve(m)
    
    final = inv %*% rep(1, dim(m)[1])
    final
}



    
