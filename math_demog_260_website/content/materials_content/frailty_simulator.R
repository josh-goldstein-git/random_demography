## A simulation of mortality frailty
## (1) get the z's: specify and draw from frailty distribution
## (2) get the ages at death from h0*z
## (3) show pop survival, hazard, and average z.
## (4) summary graphic

#######################
## (1) get the z's:  ##
#######################
## specify and draw from frailty distribution

million = 10^6
N <-  million

## uniform
## w <- .3 ## try smaller if you want
## z <- runif(N, min = 1 - w , max = 1 + w)
## mean(z)
## sd(z)
## hist(z)
## Here's code for some gamma and U-shaped beta
## gamma
my.sd <- .5
 sigma.sq <- my.sd^2
  z <- rgamma(N, shape = 1/sigma.sq, scale = sigma.sq)
  mean(z)
  sd(z)
  hist(z)
## ## beta (U-shaped)
##  z <- rbeta(N, shape1 = .5, shape2 = .5)
##  hist(z)
 ## other choices?

#########################################
## (2) simulate ages at death from h0*z ##
#########################################

## note: we call the continuous ages of death "y"
## but we'll make a table of deaths at age "x" and
## using life table notation call the count "Dx"

## we're going to use the gompertz as our baseline

source("~/Documents/teaching/math_demog_dem260/gomp_funs.R")

base.b <- 1/9
base.a <- 10^-4
y <- rgomp(N,
           b = base.b, ## doesn't vary
           a = base.a * z) ## multiplicative fixed frailty
hist(y)

###################################################
## (3) show pop survival, hazard, and average z. ##
###################################################


## first define age at death as floor(y) and then
## make a table of deaths at each age ("Dx")
get.Dx <- function(y)
{
    ## counts number of x in single year age groups
    ## including zeros when there's no one
    ## (note: built-in "table()" won't do this :(
    x <- 0:max(floor(y))
    y.fac <- factor(floor(y), levels = x)
    Dx <- tabulate(y.fac)
    names(Dx) <- x
    return(Dx)
}

Dx <- get.Dx(y)
x <- as.numeric(names(Dx))
## get lx by reverse-survival
lx <- rev(cumsum(rev(Dx)))
## get person-years as average of adjacent lx
lxpn <- c(lx[-1], 0)
Lx <- (lx + lxpn)/2
## get hazards
mx <- Dx/Lx


## (4) summary graphic

lx.base <- N * (1- pgomp(x, b = base.b, a = base.a))
Dx.base <- round(-diff(c(lx.base,0)))
mx.base <- base.a * exp(base.b * (x + .5)) ## x + .5

par(mfrow = c(2,2))
## frailty
hist(z, main = "frailty distribution")
## deaths
plot(x, Dx, main = "death distribution",
     type = "l", ylim = c(0, max(Dx.base)))
lines(x, Dx.base, lty = 2)
legend("topleft", legend = c("pop", "baseline"),
       lty = c(1, 2), cex = .8,  bty = "n")
## survival curve
plot(x, lx, main = "survival curve",
     type = "l")
lines(x, lx.base, lty = 2)
legend("topright", legend = c("pop", "baseline"),
       lty = c(1, 2), cex = .8, bty = "n")
## hazard curve
plot(x, mx, main = "hazards",
     type = "l", log = "y",
     ylim = range(mx.base, na.rm = T, finite = TRUE))
lines(x, mx.base, lty = 2)
legend("topleft", legend = c("pop", "baseline"),
       lty = c(1, 2), cex = .8, bty = "n")






