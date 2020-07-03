## "animation" of asfr change from about 2000 to about 2018

## age specific fertility rates by birth order for all countries and times
## RR means "rectangles" on Lexis surface
dt <- fread("~/Documents/hfd/Files/zip_w/asfrRRbo.txt")
dt <- dt[Code == "USA"] ## keep only US
## keep only ages 15 to 49
dt <- dt[Age %in% 15:44]
print(dt)
## put all order fertiility into a matrix
fat <- dt[, xtabs(ASFR ~ Age + Year)]
fat <- as.matrix(unclass(fat)) ## may cause problems, we'll see
fat1 <- dt[, xtabs(ASFR1 ~ Age + Year)]
fat1 <- as.matrix(unclass(fat1)) ## may cause problems, we'll see


year.vec <- colnames(fat)
age.vec <- rownames(fat)
my.year.vec <- 1975:2017
pdf("fat_movie.pdf")
for (i in 1:length(my.year.vec))
{
    my.year <- my.year.vec[i]
    plot(age.vec, fat[, paste(my.year)],
         ylim = c(0, .15),
         xlim = c(15,49),
         xlab = "age",
         ylab = "age-specific fertility rate",
         type = "l")
    title(my.year)
}
dev.off()


## now fit mixing model and redo the animation

library(data.table)
library(mixtools)
source("~/Documents/admin/socsec/analysis/tempo_mixed_functions.R")
source("~/Documents/tempo_repos/optim_clean/tempo_functions.R")

## takes few minutes to run
my.fat <- fat[, paste(my.year.vec)]
out <- get.mixed.tfr.star(my.fat)
##
out.all <- out

if (0) {
mu.mat <- get.coefs.mixed(out.all$fert.fit.list.variable.sigma)$mu.mat
lambda.mat <- get.coefs.mixed(out.all$fert.fit.list.variable.sigma)$lambda.mat
sigma.mat <- get.coefs.mixed(out.all$fert.fit.list.variable.sigma)$sigma.mat
}
mu.mat <- get.coefs.mixed(out.all$fert.fit.list)$mu.mat
lambda.mat <- get.coefs.mixed(out.all$fert.fit.list)$lambda.mat
sigma.mat <- get.coefs.mixed(out.all$fert.fit.list)$sigma.ma

matplot(my.year.vec, t(mu.mat))
abline(v = 2015)## problem here
points(c(2015, 2015), c(21.5, 30.3))

matplot(my.year.vec, t(lambda.mat))
abline(v = 2015)## problem here
points(c(2015, 2015), c(21.5, 30.3))

## interpolate 1915

colnames(lambda.mat) <- my.year.vec
colnames(mu.mat) <- my.year.vec
lambda.mat[,"2015"] <- (lambda.mat[,"2014"] + lambda.mat[,"2016"])/2
mu.mat[,"2015"]     <- (mu.mat[,"2014"]     + mu.mat[,"2016"])/2



## ## a bit faster
## my.fat.1 <- fat1[, paste(my.year.vec)]
## out1 <- get.mixed.tfr.star(my.fat.1)
## ##
## out.1 <- out1
## ## mu.mat <- get.coefs.mixed(out.1$fert.fit.list.variable.sigma)$mu.mat
## ## lambda.mat <- get.coefs.mixed(out.1$fert.fit.list.variable.sigma)$lambda.mat
## ## sigma.mat <- get.coefs.mixed(out.1$fert.fit.list.variable.sigma)$sigma.mat

## mu.mat <- get.coefs.mixed(out.1$fert.fit.list)$mu.mat
##  lambda.mat <- get.coefs.mixed(out.1$fert.fit.list)$lambda.mat
##  sigma.mat <- get.coefs.mixed(out.1$fert.fit.list)$sigma.m


## ## Aall <- my.fat.1


## plot normalized plot with the normal distributions inside.
plot.fun3 <- function(last.year, A = Aall)
{
    year.vec <- colnames(A)
    fx <- A[,paste(last.year)]
    fx <- fx/sum(fx)
    x <- as.numeric(names(fx))
    ## par(mfrow = c(1,1))
    plot(x, fx,
         ylim = c(0, .1),
         ylab = "normalized fx")
    s <- year.vec == last.year
    print(s)
    dim(mu.mat)
    length(s)
    print(x)
    print(mu.mat[1,s])
    fx1.hat <- dnorm(x, mean = mu.mat[1,s], sd = sigma.mat[1,s]) *
        lambda.mat[1,s]
    lines(x, fx1.hat, col = "red")
    fx2.hat <- dnorm(x, mean = mu.mat[2,s], sd = sigma.mat[2,s]) *
        lambda.mat[2,s]
    lines(x, fx2.hat, col = "blue")
    lines(x, fx1.hat + fx2.hat)
    title(last.year)
}


plot.fun.nonorm.3 <- function(last.year, A = Aall)
{
    ## don't normalize, but this means we have to multiply mixture by TFR
    year.vec <- colnames(A)
    fx <- A[,paste(last.year)]
##     fx <- fx/sum(fx)
    x <- as.numeric(names(fx))
    ## par(mfrow = c(1,1))
    plot(x, fx,
         ylim = c(0, .3),
         ylab = "normalized fx")
    s <- year.vec == last.year
    this.tfr <- sum(fx)
    fx1.hat <- dnorm(x, mean = mu.mat[1,s], sd = sigma.mat[1,s]) *
        lambda.mat[1,s]* this.tfr
    lines(x, fx1.hat, col = "red")
    fx2.hat <- dnorm(x, mean = mu.mat[2,s], sd = sigma.mat[2,s]) *
        lambda.mat[2,s] * this.tfr
    lines(x, fx2.hat, col = "blue")
    lines(x, fx1.hat + fx2.hat)
    title(last.year)
}

Aall <- my.fat
pdf("fat_mix_movie.pdf")
for (i in 1:length(my.year.vec))
{
    my.year <- my.year.vec[i]
    plot.fun.nonorm.3(my.year)
}
dev.off()


year.vec <- my.year.vec
A1 <- fat
##pdf("fat_1_mix_movie.pdf")
for (i in 1:length(my.year.vec))
{
    my.year <- my.year.vec[i]
    print(my.year)
    plot.fun.nonorm.3(my.year, A = A1)
}
## dev.off()


## let's do tempo adjustment

rt.mat <- t(apply(mu.mat, 1, center.diff))
tfr.vec <- apply(my.fat, 2, sum)
tfr.mat <- lambda.mat * tfr.vec
par(mfrow = c(1,2))
matplot(my.year.vec, t(tfr.mat), ylim = c(0, 3))
tfr.star.mat <- tfr.mat / (1 - rt.mat)
matplot(my.year.vec, t(tfr.star.mat), ylim = c(0,3))

tfr.star.vec <- colSums(tfr.star.mat)

par(mfrow = c(1,1))
plot(my.year.vec, tfr.vec, type = "l",
     ylim = c(1, 3))

lines(my.year.vec, tfr.star.vec, lty = 2)

