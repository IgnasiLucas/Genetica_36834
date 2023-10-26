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
         warning('Interpreting input as genotypes from one or more loci.')
      } else {
         stop('Wrong input. Use a list of length-two vectors for genotypes in one or more loci.')
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
   # gametes are character vectors, with genes assumed in same order!
   stopifnot(length(gamete1) == length(gamete2))
   stopifnot(is.character(gamete1))
   stopifnot(is.character(gamete2))
   genotypes <- mapply(c, gamete1, gamete2, SIMPLIFY = FALSE)
   zigot <- diploid(genotypes)
   return(zigot)
}

is.diploid <- function(x) {
   class(x) == 'diploid'
}

'*.diploid' <- function(mum, dad) {
   stopifnot(is.diploid(mum) && is.diploid(dad))
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

phenotype <- function(x, map = 'dihybrid'){
   # Assumes every gene uses one different letter, with
   # upper case meaning "dominant" and lower case, "recessive".
   if (is.matrix(x)) {
      # In case x is the output of cross(), main intended use.
      z <- suppressWarnings(lapply(strsplit(x, '[, ]'), diploid))
   } else if (all(sapply(x, class) == 'diploid')) {
      z <- x
   } else if (class(x) == 'diploid') {
      z <- list(x)
   } else stop('Wrong input.')
   # Translates genotypes to numbers: 
   f <- sapply(z, function(y) {
      sum(sapply(y, function(w) length(grep(paste(w, collapse = '|'), LETTERS))) * 2^(0:(length(y)-1)))
   })
   map <- match.arg(map, c('dihybrid','simple.recessive','simple.dominant', 'quantitative',
                           'double.recessive','double.dominant','double.dominant.recessive'))
   if (map == 'dihybrid') {
      f <- f
   } else if (map == 'simple.recessive') {
      f[f == 2] <- 0
   } else if (map == 'simple.dominant') {
      f[f == 1] <- 3
   } else if (map == 'double.recessive') {
      f[f == 1] <- 0
      f[f == 2] <- 0
   } else if (map == 'double.dominant') {
      f[f == 1] <- 3
      f[f == 2] <- 3
   } else if (map == 'double.dominant.recessive') {
      f[f == 1] <- 3
      f[f == 0] <- 3
   } else if (map == 'quantitative') {
      f <- sapply(z, function(y) {
         sum(sapply(y, function(w) sum(w %in% LETTERS)))  
      })
   }
   return(structure(f, dim = dim(x)))
}

punnet <- function(x, palette = 'Egypt', map = 'dihybrid') {
   stopifnot(is.matrix(x))
   stopifnot(is.character(x))
   Map <- match.arg(map, c('dihybrid', 'simple.recessive', 'simple.dominant', 'quantitative',
                           'double.recessive', 'double.dominant', 'double.dominant.recessive'))
   n <- dim(x)[1]
   m <- dim(x)[2]
   Genotips <- x
   Fenotips <- phenotype(x, map = Map) + 1
   par(mar = c(1, 6, 6, 1))
   image(x = 1:m, y = 1:n, z = t(Fenotips[n:1,]),
         col = MetBrewer::met.brewer(palette, max(Fenotips), direction = 1),
         axes = FALSE, xlab = '', ylab = '')
   axis(2, at = seq(from = 0.5, to = 0.5 + n, by = 1),
        labels = FALSE)
   axis(2, at = 1:n, tick = FALSE,
        labels = row.names(Genotips)[n:1], cex.axis = 2, las = 2)
   axis(3, at = seq(from = 0.5, to = 0.5 + m, by = 1),
        labels = FALSE)
   axis(3, at = 1:m, tick = FALSE,
        labels = colnames(Genotips), cex.axis = 2)
   abline(h = seq(from = 0.5, to = 0.5 + m, by = 1))
   abline(v = seq(from = 0.5, to = 0.5 + n, by = 1))
   text(x = rep(1:m, each = n), y = rep(n:1, m), labels = Genotips, cex = 32/(m * m))
   mtext('First gamete', side = 2, line = 4, cex = 2)
   mtext('Second gamete', side = 3, line = 4, cex = 2)
   par(mar = c(5, 4, 4, 2))
}

histo <- function(x, map = 'dihybrid', palette = 'Egypt') {
   stopifnot(is.matrix(x))
   stopifnot(is.character(x))
   Map <- match.arg(map, c('dihybrid', 'simple.recessive', 'simple.dominant',
                           'double.recessive', 'double.dominant', 'double.dominant.recessive'))
   fenotips <- suppressWarnings(phenotype(x, map = Map)) + 1
   F <- as.vector(table(fenotips))
   barplot(F, col = MetBrewer::met.brewer(palette, max(fenotips), direction = 1),
           xlab = 'Phenotypes', ylab = 'Expected frequency',
           main = paste('Proportions', paste0(sort(F, decreasing = TRUE), collapse = ':')))
}
