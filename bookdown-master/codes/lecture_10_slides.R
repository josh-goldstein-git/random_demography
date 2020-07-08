## -----------------------------------------------------------------------------
fwm <- function(N, n_gen, mu = 0) ## mu != 4/N
{
    ## simulate fisher-wright (with mutations)
    x <- paste(1:N) ## starting types
    A <- matrix(NA, nrow = n_gen, ncol = N)
    for (i in 1:n_gen)
    {
        A[i,] <- x
        x <- sample(x, size = N, replace = T)
        x <- mut(x, mu)
        x
    }
    return(A) ## matrix of types, each line a generation.
}


## -----------------------------------------------------------------------------
mut <- function(x, mu)
{
    ## m, the individuals that mutate
    m <- which(rbinom(length(x), 1, mu) == 1)
    if (length(m) == 0) ## if no-one mutates
        return(x)
    ## add a suffix to their ID, so it will be unique (infinite alleles)
    suffix <- 10000*round(runif(length(m)),4)
    x[m] <- paste0(x[m], ".", suffix)
    x
}



## -----------------------------------------------------------------------------
set.seed(1)
fwm(N = 10, n_gen = 2, mu = 0)
fwm(N = 10, n_gen = 2, mu = 0)


## ----size="tiny", fig.height = 3.5--------------------------------------------
set.seed(1)
A <- fwm(N = 10, n_gen = 20, mu = 0)
tt <- table(A, row(A)) ## count types by row
ptt <- prop.table(tt, 2) ## proportions
matplot(t(ptt), type = 'l', lty = 1, main = "FW simu")
text(x = 4, y = jitter(ptt[,4]), rownames(ptt), col = 1:6)


## ----size = "tiny", fig.height = 4--------------------------------------------
set.seed(1)
A <- fwm(N = 100, n_gen = 200, mu = 0)
tt <- table(A, row(A)) ## count types by row
ptt <- prop.table(tt, 2) ## proportions
matplot(t(ptt), type = 'l', lty = 1)


## ----echo=F, fig.height = 5---------------------------------------------------
set.seed(1)
A <- fwm(N = 10, n_gen = 20, mu = 0)
tt <- table(A, row(A)) ## count types by row
ptt <- prop.table(tt, 2) ## proportions
matplot(t(ptt), type = 'l', lty = 1, main = "FW simu")
text(x = 4, y = jitter(ptt[,4]), rownames(ptt), col = 1:6)


## -----------------------------------------------------------------------------
T.vec <- NULL
all.the.same <- function(x){length(unique(x)) == 1}
set.seed(10)
for (i in 1:100)
{
    A <- fwm(N = 100, n_gen = 1000,mu = 0)
    extinction_time = min(row(A)[apply(A, 1, all.the.same)])
    T.vec[i] <- extinction_time
}
mean(T.vec)

