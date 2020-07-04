## get center.diff

source("/hdir/0/fmenares/Book/bookdown-master/codes/utility_functions.R")

## fit mixed model to series of periods

fit.mixed.from.A <- function(A, disn = "normal")
{
    n_years = ncol(A)
    year.vec <- colnames(A)
    age <- as.numeric(rownames(A))
    fert.fit.list <- vector("list", length = n_years)
    for (i in 1:n_years)
    {
        this.year <- year.vec[i]
        print(this.year)
        fx <- A[,this.year]
        z <- rep(age, round(fx * 10^4))
        if (disn == "normal")
            fert.fit <- normalmixEM(z, lambda = c(.2, .8),
                                    mu = c(20, 30),
##                                    fast = TRUE,
                                    maxit = 5000,
                                    sigma = c(3, 5))
        if (disn == "gamma")
            fert.fit <- gammamixEM(z, lambda = c(.2, .8),
                                   alpha = c(60, 30),
                                   beta = c(.3, .8))
        fert.fit.list[[i]] <- fert.fit
    }
    return(fert.fit.list)
}

## extract coefs from fert.fit.list

get.coefs.mixed <- function(fert.fit.list, disn = "normal")
{
    L = fert.fit.list
    if (disn == "normal") {
    mu.mat <- sapply(L, "[[", "mu")
    sigma.mat <- sapply(L, "[[", "sigma")
    lambda.mat <- sapply(L, "[[", "lambda")
    mixed.coefs.list <-
    list(mu.mat = mu.mat,
         sigma.mat = sigma.mat,
         lambda.mat = lambda.mat)
    }
    if (disn == "gamma")
    {
    alpha.mat <- sapply(L, "[[", "alpha")
    beta.mat <- sapply(L, "[[", "beta")
    lambda.mat <- sapply(L, "[[", "lambda")
    mixed.coefs.list <-
        list(alpha.mat = alpha.mat,
             beta.mat = beta.mat,
             lambda.mat = lambda.mat)
        }
    return(mixed.coefs.list)
}

fit.mixed.again.constant.sigma <- function(A, fert.fit.list,
                                           constant.sigma.vec = NULL)
{
    ## set sigma at average
    if (is.null(constant.sigma.vec))
        {
            sigma.mat <- get.coefs.mixed(fert.fit.list)$sigma.mat
            constant.sigma.vec <- apply(sigma.mat, 1, mean)
            ## 1 sigma for each component
        }

    ## now re-estimate using old estimats as starting values
    ## along with constant sigma
    n_years = ncol(A)
    year.vec <- colnames(A)
    age <- as.numeric(rownames(A))
    fert.fit.list.constant.sigma <- vector("list", length = n_years)
    for (i in 1:n_years)
    {
        this.year <- year.vec[i]
        print(this.year)
        fx <- A[,this.year]
        z <- rep(age, round(fx * 10^4))
        est <- fert.fit.list[[i]]
        fert.fit <- normalmixEM(z, lambda = est$lambda,
                                mu = est$mu,
##                                fast = TRUE,
                                maxit = 5000,
                                sd.constr  = constant.sigma.vec)
        fert.fit.list.constant.sigma[[i]] <- fert.fit
    }
    return(fert.fit.list.constant.sigma)
}

## visualize fits

show.fits <- function(A, fert.fit.list)
{

    return(NULL)
}

## get tempo-adjusted fertility

get.tfr.star <- function(A, fert.fit.list)
{
    mu.mat <- get.coefs.mixed(fert.fit.list)$mu.mat
    lambda.mat <- get.coefs.mixed(fert.fit.list)$lambda.mat

    rt.a <- center.diff(mu.mat[1,], end.fill = T)
    rt.b <- center.diff(mu.mat[2,], end.fill = T)
    tfr <- apply(A, 2, sum)
    tfr.a <- lambda.mat[1,] * tfr
    tfr.b <- lambda.mat[2,] * tfr
    tfr.a.star <- tfr.a/(1 - rt.a)
    tfr.b.star <- tfr.b/(1 - rt.b)
    tfr.ab.star <- tfr.a.star + tfr.b.star
    tfr.star <- tfr.ab.star
    return(tfr.star)
}

## the whole shebang

get.mixed.tfr.star <- function(A)
{
    fert.fit.list <- fit.mixed.from.A(A)
    fert.fit.list.constant.sigma <-
        fit.mixed.again.constant.sigma(A, fert.fit.list)
    mixed.coefs.list <- get.coefs.mixed(fert.fit.list)
    tfr.star.variable.sigma <- get.tfr.star(A, fert.fit.list)
    tfr.star.constant.sigma <-
        get.tfr.star(A, fert.fit.list.constant.sigma)
    out <- list(
        tfr.star.variable.sigma = tfr.star.variable.sigma,
        tfr.star = tfr.star.constant.sigma,
        fert.fit.list.variable.sigma = fert.fit.list,
        fert.fit.list = fert.fit.list.constant.sigma)
    return(out)
}

## let's try various sigmas




if(0) {

## read in data
dt <- fread("~/Documents/hfd/Files/zip_w/asfrRRbo.txt")
dt[, ASFR3p := ASFR3 + ASFR4 + ASFR5p]
dt[Age == "12-", Age := 12]
dt[Age == "55+", Age := 55]
dt[, x := as.numeric(Age)]
dt[, x := as.numeric(x)]
dt <- dt[Code == "USA" & Year >= 2005]


A1 <- xtabs(ASFR1 ~ x + Year, data = dt) ## parity 1
A <- A1
fert.fit.list <- fit.mixed.from.A(A)
fert.fit.list.constant.sigma <- fit.mixed.again.constant.sigma(A, fert.fit.list)
mixed.coefs.list <- get.coefs.mixed(fert.fit.list)
tfr.star <- get.tfr.star(A, fert.fit.list)
tfr.star.constant.sigma <- get.tfr.star(A, fert.fit.list.constant.sigma)
tfr.star.1 <- tfr.star.constant.sigma


}
