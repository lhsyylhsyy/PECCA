setwd('C:\\Users\\Ning\\Dropbox\\Research\\SIR\\ECCA simulation 08262018')
# ECCA example 08262018
n <- 250
set.seed(7202016)
x12 <- matrix(rnorm(2*n, 0, 1), n, 2)
e1 <- rnorm(n, 0, 1)
e2 <- rnorm(n, 0, 1)
y23 <- matrix(rnorm(2*n, 0, 1), n, 2)

x3 <- 0.8 * log(abs(y23[,1] + y23[,2]))  + 0.4 * e2
y1 <- 0.3 * (x12[,1] + x12[,2])^2 + 0.2 * e1

x <- cbind(x12, x3)
y <- cbind(y1, y23)

model1 <- ECCA(x, y, 2)
model.cca <- cc(x, y)

x.variates <- x %*% model1[[2]]
y.variates <- y %*% model1[[1]]

x.variates.cca <- x %*% model.cca$xcoef
y.variates.cca <- y %*% model.cca$ycoef

par(mfrow=c(2,2))
plot(x.variates[,1], y.variates[,1],
     xlab = 'X variate 1', ylab = 'Y variate 1', main = '(a)')
plot(x.variates[,2], y.variates[,2],
     xlab = 'X variate 2', ylab = 'Y variate 2', main = '(b)')
plot(x.variates.cca[,1], y.variates.cca[,1],
     xlab = 'X variate 1', ylab = 'Y variate 1', main = '(c)')
plot(x.variates.cca[,2], y.variates.cca[,2],
     xlab = 'X variate 2', ylab = 'Y variate 2', main = '(d)')


