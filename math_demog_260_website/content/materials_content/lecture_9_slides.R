## ----echo = T, size = "small"-------------------------------------------------
branch <- function(n_max = 30, p_k = c(p0, p1, p2), Z1 = 1)
{
    ## note: this returns 0s when extinct
    k <- 0:(length(p_k)-1)
    Z.vec <- rep(NA, n_max)
    Z.vec[1] <- Z1
    for (i in 1:(n_max-1))
    {
        Z.vec[i+1] <- sum(sample(x = k,
                                 size = Z.vec[i],
                                 replace = T,
                                 prob = p_k))
    }
    return(Z.vec)
}
p0 = .3; p1 = .4; p2 = .3 ## what is m?


## ----echo = T, eval = F-------------------------------------------------------
## n_trials = 1000; n_gen = 100
## Z.mat <- matrix(NA, n_trials, n_gen)
## set.seed(131)
## for (i in 1:n_trials)
##     Z.mat[i,] <- branch(n_max = n_gen)
## matplot(t(Z.mat),
##         type = "l", lty = 1, ylab = "Zn", xlab = "n")


## ----echo = F, eval = T, fig.height = 5---------------------------------------
n_trials = 1000; n_gen = 100
Z.mat <- matrix(NA, n_trials, n_gen)
set.seed(131)
for (i in 1:n_trials)
    Z.mat[i,] <- branch(n_max = n_gen)
matplot(t(Z.mat),
        type = "l", lty = 1, ylab = "Zn", xlab = "n")


## ----echo = T, eval = F-------------------------------------------------------
## Zn_bar = apply(Z.mat, 2, mean)
## n <- 1:ncol(Z.mat)
## proportion.zero <- function(x){prop.table(table(x == 0))["TRUE"]}
## d_n = apply(Z.mat, 2, proportion.zero)
## Z.mat.na <- Z.mat; Z.mat.na[Z.mat == 0] <- NA
## Zn_surv_bar = apply(Z.mat.na, 2, mean, na.rm = T)
## par(mfrow = c(1,3))
## plot(n, Zn_bar, main = "Mean Zn")
## plot(n, d_n, main = "Fraction extinct")
## plot(n, Zn_surv_bar)
## ## insert code here for Zn_surv_bar.hat and add a line


## ----echo = F, eval = T, fig.height = 5---------------------------------------
Zn_bar = apply(Z.mat, 2, mean)
n <- 1:ncol(Z.mat)
proportion.zero <- function(x){prop.table(table(x == 0))["TRUE"]}
d_n = apply(Z.mat, 2, proportion.zero)
Z.mat.na <- Z.mat; Z.mat.na[Z.mat == 0] <- NA
Zn_surv_bar = apply(Z.mat.na, 2, mean, na.rm = T)
par(mfrow = c(1,3))
plot(n, Zn_bar, main = "Mean Zn")
plot(n, d_n, main = "Fraction extinct")
plot(n, Zn_surv_bar)
## insert code here for Zn_surv_bar.hat and add a line


## ----fig.height = 5-----------------------------------------------------------
var_Zn = apply(Z.mat, 2, var)
n <- 1:ncol(Z.mat)
plot(n, var_Zn)


## ----fig.height = 4-----------------------------------------------------------
b = 0.2126 ;   c = 0.5893
kk = 1:10  ;   p_kk = b * c^(kk-1)
p0 = b/(1-c)
k = c(0, kk) ;   p_k = c(p0, p_kk)
plot(k, p_k)


## ----eval = F, echo = T, size = "tiny"----------------------------------------
## ## b = 0.2126 ;   c = 0.5893 ## lotka
## b = 0.3666; c = .5533 ## Keyfitz (from GS)
## m = b / (1-c)^2 ## [1] 1.260416
## d = (1 - b - c) / (c * (1-c))           #[1] 0.8185088
## par(mfrow = c(1,1))
## for (i in 1:3)
## {
##     n = i
##     p0_n = d * (m^n - 1)/ (m^n -d)
##     j = kk
##     pj_n = m^n *
##         ((1-d) / (m^n - d))^2 *
##         ((m^n - 1)/(m^n - d))^(j-1)
##     pk_n <- c(p0_n, pj_n)
##     if (i == 1)
##         plot(k, pk_n, type = "l", log = "")
##     if (i > 1)
##         lines(k, pk_n, col = i)
## }


## ----eval = T, echo = F, fig.height = 5---------------------------------------
## b = 0.2126 ;   c = 0.5893 ## lotka
b = 0.3666; c = .5533 ## Keyfitz (from G&S)
m = b / (1-c)^2 ## [1] 1.260416
d = (1 - b - c) / (c * (1-c))           #[1] 0.8185088
par(mfrow = c(1,1))
for (i in 1:3)
{
    n = i
    print(n)
    p0_n = d * (m^n - 1)/ (m^n -d)
    j = kk
    pj_n = m^n * ((1-d) / (m^n - d))^2 * ((m^n - 1)/(m^n - d))^(j-1)
    pk_n <- c(p0_n, pj_n)
    if (i == 1)
        plot(k, pk_n, type = "l", log = "y")
    if (i > 1)
        lines(k, pk_n, col = i)
}

