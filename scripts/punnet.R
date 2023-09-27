monohibrid <- function(n=0, palette='Egypt'){
   punnet <- matrix(c(1, 0, 1, 1), ncol = 2, byrow = FALSE)
   par(mar = c(1, 6, 6, 1))
   image(x = c(1, 2), y = c(1, 2), z = punnet,
         col = MetBrewer::met.brewer(palette, 2),
         axes = FALSE, xlab = '', ylab = '')
   axis(2, at = c(0.5, 1.5, 2.5), labels = c('', '', ''))
   axis(2, at = c(1, 2), tick = FALSE, labels = c('a', 'A'), las = 2, cex.axis = 2)
   axis(3, at = c(0.5, 1.5, 2.5), labels = c('', '', ''))
   axis(3, at = c(1, 2), tick = FALSE, labels = c('A', 'a'), cex.axis = 2)
   abline(h = c(0.5, 1.5, 2.5))
   abline(v = c(0.5, 1.5, 2.5))
   text(x = c(1, 1, 2, 2), y = c(1, 2, 1, 2),
        labels = c('aA', 'AA', 'aa', 'Aa'), cex = 2)
   mtext('Gameto masculino', side = 2, line = 4, cex = 2)
   mtext('Gameto femenino',  side = 3, line = 4, cex = 2)
}

dihibrid <- function(n=0, palette = 'Egypt'){
   punnet <- matrix(c(3, 2, 1, 0, 3, 3, 1, 1, 3, 2, 3, 2, 3, 3, 3, 3),
                    ncol = 4, byrow = FALSE)
   par(mar = c(1, 6, 6, 1), mfrow = c(1, 2))
   image(x = 1:4, y = 1:4, z = punnet,
         col = MetBrewer::met.brewer(palette, 4, direction = 1),
         axes = FALSE, xlab = '', ylab = '')
   axis(2, at = c(0.5, 1.5, 2.5, 3.5, 4.5), labels = c('', '', '', '', ''))
   axis(2, at = c(1, 2, 3, 4), tick = FALSE,
        labels = c('ab', 'aB', 'Ab', 'AB'), cex.axis = 2, las = 2)
   axis(3, at = c(0.5, 1.5, 2.5, 3.5, 4.5), labels = c('', '', '', '', ''))
   axis(3, at = c(1, 2, 3, 4), tick = FALSE,
        labels = c('AB', 'Ab', 'aB', 'ab'), cex.axis = 2)
   abline(h = c(0.5, 1.5, 2.5, 3.5, 4.5))
   abline(v = c(0.5, 1.5, 2.5, 3.5, 4.5))
   text(x = c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4),
        y = c(1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4),
        labels = c('aA bB', 'aA BB', 'AA bB', 'AA BB',
                   'aA bb', 'aA Bb', 'AA bb', 'AA Bb',
                   'aa bB', 'aa BB', 'Aa bB', 'Aa BB',
                   'aa bb', 'aa Bb', 'Aa bb', 'Aa Bb'),
        cex = 2)
   mtext('Gameto masculino', side = 2, line = 4, cex = 2)
   mtext('Gameto femenino',  side = 3, line = 4, cex = 2)

   par(mar = c(5, 4, 4, 2))
   fenotips <- as.vector(table(punnet))
   fenotips <- sort(fenotips, decreasing = TRUE)
   barplot(fenotips, col = MetBrewer::met.brewer(palette, 4, direction = -1),
           names.arg = c('A-B-', 'A-bb', 'aaB-', 'aabb'),
           xlab = 'Clases fenotípicas', ylab = 'Frecuencia esperada',
           main = paste('Proporciones', paste0(fenotips, collapse = ':')))
}

trihibrid <- function(n=0, palette = 'Egypt'){
   punnet <- matrix(c(8, 7, 6, 5, 4, 3, 2, 1, 8, 8, 6, 5, 6, 5, 2, 2,
                      8, 7, 8, 5, 7, 3, 5, 3, 8, 7, 6, 8, 4, 7, 6, 4,
                      8, 8, 8, 5, 8, 5, 5, 5, 8, 8, 6, 8, 6, 8, 6, 6,
                      8, 7, 8, 8, 7, 7, 8, 7, 8, 8, 8, 8, 8, 8, 8, 8),
                    ncol = 8, byrow = FALSE)
   par(mar = c(1, 6, 6, 1), mfrow = c(1, 2))
   image(x = 1:8, y = 1:8, z = punnet,
         col = MetBrewer::met.brewer(palette, 8, direction = 1),
         axes = FALSE, xlab = '', ylab = '')
   axis(2, at = c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5), labels = rep('', 9))
   axis(2, at = c(1, 2, 3, 4, 5, 6, 7, 8), tick = FALSE,
        labels = c('abc', 'abC', 'aBc', 'Abc', 'aBC', 'AbC', 'ABc', 'ABC'), cex.axis = 1, las = 2)
   axis(3, at = c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5), labels = rep('', 9))
   axis(3, at = c(1, 2, 3, 4, 5, 6, 7, 8), tick = FALSE,
        labels = c('ABC', 'ABc', 'AbC', 'aBC', 'Abc', 'aBc', 'abC', 'abc'), cex.axis = 1)
   abline(h = c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5))
   abline(v = c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5))
#   text(x = c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4),
#        y = c(1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4),
#        labels = c('aA bB', 'aA BB', 'AA bB', 'AA BB',
#                   'aA bb', 'aA Bb', 'AA bb', 'AA Bb',
#                   'aa bB', 'aa BB', 'Aa bB', 'Aa BB',
#                   'aa bb', 'aa Bb', 'Aa bb', 'Aa Bb'),
#        cex = 2)
   mtext('Gameto masculino', side = 2, line = 4, cex = 2)
   mtext('Gameto femenino',  side = 3, line = 4, cex = 2)

   par(mar = c(6, 4, 4, 2))
   fenotips <- as.vector(table(punnet))
   fenotips <- sort(fenotips, decreasing = TRUE)
   barplot(fenotips, col = MetBrewer::met.brewer(palette, 8, direction = -1),
           names.arg = c('A-B-C-', 'A-B-cc', 'A-bbC-', 'aaB-C-',
                         'A-bbcc', 'aaB-cc', 'aabbC-', 'aabbcc'),
           las = 2,
           ylab = 'Frecuencia esperada',
           main = paste('Proporciones', paste0(fenotips, collapse = ':')))
   mtext('Clases fenotípicas', side = 1, line = 4.5)
}
