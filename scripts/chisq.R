library(MetBrewer)
x <- seq(from = 0, to = 8, by = 0.01)
y <- sapply(1:5, function(z) dchisq(x, z))
p <- met.brewer('Egypt', 5)
plot(c(0,8), c(0, 0.5), type = 'n', xlab = 'X', ylab = '',
     main = expression(paste('Densidad de probabilidad ', chi[k]^2)))
lapply(1:5, function(z) lines(x, y[,z], col = p[z], lwd = 2))
legend(6, 0.45, paste('k', 1:5, sep = '='), col = p, lwd = 2)
