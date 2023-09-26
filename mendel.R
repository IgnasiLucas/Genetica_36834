new_diploid <- function(x, chr = 1:length(x), pos = rep(1, length(x))){
   stopifnot(all(sapply(x, is.character)))
   stopifnot(is.list(x))
   stopifnot(all(sapply(x, length) == 2))
   structure(x, chr = chr, pos = pos, class = 'diploid')
}
diploid <- function(x, chr = 1:length(x), pos = rep(1, length(x))){
   if (is.character(x)) {
      if (length(x) == 2) {
         if (all(sapply(x, nchar) == 1)) {
            x <- list(l1 = x)
            warning('Interpreting input as two alleles from one locus.')
         } else if (all(sapply(x, nchar) == 2)) {
            x <- sapply(x, strsplit, '')
            warning('Interpreting input as genotypes from two loci.')
         }
      } else if (all(sapply(x, nchar) == 2)) {
         x <- sapply(x, strsplit, '')
         warning('Interpreting input as genotypes from >2 loci.')
      } else {
         error('Wrong input. Use a list of length-two vectors for genotypes in one or more loci.')
      }
   } else if (is.list(x) && all(sapply(x, length) == 2)) {
      x <- lapply(x, as.character)
   } else {
      error('Wrong input. Use a list length-two vectors for genotypes in one or more loci.')
   }
   names(x) <- paste0('l', 1:length(x))
   new_diploid(x)
}

meiosis <- function(x) {
   UseMethod('meiosis')
}

meiosis.diploid <- function(x) {
   stopifnot(class(x) == 'diploid')
   NumGenes <- length(x)
   NumAlleles <- sapply(lapply(x, unique), length)
   NumDifGametes <- prod(NumAlleles)
   gametes <- matrix(NA, ncol = NumGenes, nrow = NumDifGametes)
   for (i in 1:NumGenes) {
      Each <- NumDifGametes / prod(NumAlleles[1:i])
      gametes[, i] <- rep(sort(unique(x[[i]]), decreasing = TRUE), each = Each)
   }
   gametes <- lapply(1:NumDifGametes, function(z) gametes[z, ])
   return(gametes)
}

fusion <- function(gamete1, gamete2) {
   stopifnot(length(gamete1) == length(gamete2))
   stopifnot(is.character(gamete1))
   stopifnot(is.character(gamete2))
   genotypes <- mapply(c, gamete1, gamete2, SIMPLIFY = FALSE)
   zigot <- diploid(genotypes)
   return(zigot)
}

cross <- function(mum, dad) {
   common.genes   <- intersect(names(mum), names(dad))
   stopifnot(length(common.genes) > 0)
   mum <- diploid(mum[common.genes])
   dad <- diploid(dad[common.genes])
   female.gametes <- meiosis(mum)
   male.gametes   <- meiosis(dad)
   punnet <- matrix(0, nrow = length(female.gametes),
                    ncol = length(male.gametes))
   colnames(punnet)  <- sapply(male.gametes, paste, collapse = '')
   row.names(punnet) <- sapply(female.gametes, paste, collapse = '')
   for (g1 in female.gametes) {
      for (g2 in male.gametes) {
         zigot <- fusion(g1, g2)
         i <- paste(g1, collapse = '')
         j <- paste(g2, collapse = '')
         punnet[i, j] <- paste(
            lapply(zigot, function(x) paste(sort(x, decreasing = TRUE), collapse = '')),
            collapse = ',')
      }
   }
   return(punnet)
}

fenotip <- function(x, map = 'dihibrid'){
   # Assumes every gene uses one different letter, with
   # upper case meaning "dominant" and lower case, "recessive".
   if (is.matrix(x)) {
      # In case x is the output of cross(), main intended use.
      z <- lapply(strsplit(x, '[, ]'), diploid)
   } else if (all(sapply(x, class) == 'diploid')) {
      z <- x
   } else if (class(x) == 'diploid') {
      z <- list(x)
   } else stop('Wrong input.')
   f <- sapply(z, function(y) {
      sum(sapply(y, function(w) length(grep(paste(w, collapse = '|'), LETTERS))) * 2^(0:(length(y)-1)))
   })
   if (map == 'dihibrid') {
      return(structure(f, dim = dim(x)))
   } else if (mapp == 'simple.recesiva') {

   }
}
