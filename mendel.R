ploidy <- function(genotype){
   stopifnot(all(sapply(genotype, length) == length(genotype[[1]])))
   return(length(genotype[[1]]))   
}

meiosis <- function(genotype, n = 0){
   # genotype must be a list of genotypes of one or more genes.
   stopifnot(is.list(genotype))
   p <- ploidy(genotype)
   stopifnot(p %% 2 == 0)
   NumGenes <- length(genotype)
   NumAlleles <- sapply(lapply(genotype, unique), length)
   NumDifGametes <- prod(NumAlleles)
   gametes <- matrix(NA, ncol = NumGenes, nrow = NumDifGametes)
   for (i in 1:NumGenes) {
      Each <- NumDifGametes / prod(NumAlleles[1:i])
      gametes[, i] <- rep(sort(unique(genotype[[i]]), decreasing = TRUE), each = Each)
   }
   if (n == 0) {
      apply(gametes, 1, paste, collapse = '')
   } else {
      sample(apply(gametes, 1, paste, collapse = ''), n, replace = TRUE)
   }
}
fusion <- function(gamete1, gamete2) {
   stopifnot(toupper(gamete1) == toupper(gamete2))
   zigor <- paste0(unlist(strsplit(gamete1, '')),
                   unlist(strsplit(gamete2, '')))
   zigot <- lapply(zigor, function(x){
      paste0(sort(unlist(strsplit(x, '')), decreasing = TRUE), collapse = '')
   })
   names(zigot) <- toupper(unlist(strsplit(gamete1, '')))
   return(zigot)
}
cross <- function(mum, dad, n = 0) {
   common.genes   <- intersect(names(mum), names(dad))
   stopifnot(length(common.genes) > 0)
   mum <- mum[common.genes]
   dad <- dad[common.genes]
   female.gametes <- meiosis(mum)
   male.gametes   <- meiosis(dad)
   punnet <- matrix(0, nrow = length(female.gametes),
                    ncol = length(male.gametes))
   colnames(punnet)  <- male.gametes
   row.names(punnet) <- female.gametes
   for (g1 in female.gametes) {
      for (g2 in male.gametes) {
         punnet[g1, g2] <- paste0(fusion(g1, g2), collapse  = '')
      }
   }
   return(punnet)
}
