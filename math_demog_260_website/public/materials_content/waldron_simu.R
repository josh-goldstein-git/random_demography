## Waldron Simulation using Gamma-Gompertz

## We want to see if heterogeneity alone can explain the pattern in
## Waldron's observations of hazard ratios declining by age and
## increasing from one cohort to the next.

## The key idea is that moving from one cohort to the next, holding
## age constant, is like moving to younger ages within a cohort. Just
## as moving to younger ages makes the observed hazard ratio larger,
## moving to a lower baseline mortality rate, holding age constant,
## also makes the hazard ratio larger.


## Here is the gamma gompertz formula
## mu.bar <- a * exp(b * x) / (1 + (a * s2 / b) * (exp(b*x) - 1))
## where s2 is sigma^2

t <- 1:5 ## five different cohrots (say each 10 years apart)
k <- .25 ## amount of mortality decline per decade (about 2.5% per year)
b <- .1 ## gompertz slope
a0 <- 10^{-3} ## starting hazard (rather high)
a0.vec <- a0 * exp(-k * t)
x <- 0:100
mu.1.bar.xt <- matrix(NA, length(t), length(x))
mu.2.bar.xt <- matrix(NA, length(t), length(x))
dimnames(mu.1.bar.xt) <- list(t, x)
dimnames(mu.2.bar.xt) <- list(t, x)
s2 <- .2

for (i in 1:length(t))
{
    a1 <- a0.vec[i]
    a2 <- a1 * 2
    mu.1.bar.xt[i,] <- a1 * exp(b * x) /
        (1 + (a1 * s2 / b) * (exp(b*x) - 1))
    mu.2.bar.xt[i,] <- a2 * exp(b * x) /
        (1 + (a2 * s2 / b) * (exp(b*x) - 1))
}

R.mat <- mu.2.bar.xt/mu.1.bar.xt

waldron.simu <- R.mat[, x %in% seq(60, 100, 10)]
print(waldron.simu)
##         60       70       80       90      100
## 1 1.443725 1.226597 1.097246 1.038110 1.014365
## 2 1.505985 1.273363 1.121510 1.048411 1.018370
## 3 1.568061 1.325716 1.150818 1.061317 1.023465
## 4 1.628069 1.382813 1.185699 1.077385 1.029930
## 5 1.684373 1.443338 1.226496 1.097228 1.038107
## Interestingly, you can see that after 5 decades, the R(x) curve has shifted over almost exactly 10 years. R(90, decade 1) nearly equals R(100, decade 5) and similarly for age 60, 70, and 80.


## we now block the lower right cells (so that our table looks like Waldrons

waldron.upper <- waldron.simu
waldron.upper[2,5] <- NA
waldron.upper[3,4:5] <- NA
waldron.upper[4,3:5] <- NA
waldron.upper[5,2:5] <- NA


print(round(waldron.upper,2))
##     60   70   80   90  100
## 1 1.44 1.23 1.10 1.04 1.01
## 2 1.51 1.27 1.12 1.05   NA
## 3 1.57 1.33 1.15   NA   NA
## 4 1.63 1.38   NA   NA   NA
## 5 1.68   NA   NA   NA   NA

## We see a more rapid decline by age, and a less rapid increase by
## cohort. But qualitatively the pattern looks similar.

## Further investigation could try to choose more accurate a0 and b
## parameters to match better with observed mortality rates.

## Note: the true risk of being in group 2 instead of 1 is R = 2. So
## it's not that by cohort 5 we're seeing a distortion of the
## R(60). Rather we're seeing convergence of the population value to
## the individual risk.





