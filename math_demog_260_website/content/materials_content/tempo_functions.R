bf.fit.simple <- function(fat)
{
    ## takes matrix (e.g. 1 parity or all parity combined)
    tfr.vec <- apply(fat, 2, sum)
    ## get rt
    mu.vec <- apply(fat, 2, get.mean, x = 1:nrow(fat))
    rt.vec <- center.diff(mu.vec, end.fill = T)
    ## get bf tfr*
    tfr.star <- tfr.vec/(1-rt.vec)
    return(list(tfr.star=tfr.star, rt.vec=rt.vec))
}

bf.fit <- function(fat.array)
  {
    ## get parity specific tfr's
    tfr.mat <- apply(fat.array, 2:3, sum)
    ## get rt
    rownames(fat.array) <- remove.plusminus(rownames(fat.array))
    x <- as.numeric(rownames(fat.array))
    mu.mat <- apply(fat.array, 2:3, get.mean, x = x)
    rt.mat <- apply(mu.mat, 2, center.diff, end.fill = T)
    ## get bf tfr*
    tfr.star.mat <- tfr.mat/(1-rt.mat)
    bf.tfr.star <- apply(tfr.star.mat[,-1], 1, sum)
    tfr.star.out <- cbind(tfr.star.mat, bf.tfr.star)
    return(list(tfr.star.out=tfr.star.out, rt.mat=rt.mat))
  }


get.Sp.vec.regress <- function(fat,
                       return.se = FALSE,
                       smoothed = NULL)
{

    fac.big <- per2coh(fat)
    fac <- fac.big[,!apply(fac.big, 2, all.na)] # shrink to right size

    fa.mat <- apply(fat, 2, center.diff, end.fill = T)
    ft.mat <- t(apply(fat, 1, center.diff, end.fill = T))

    coh.names <- colnames(fac)
    coh <- coh.names[per2c(fat)]

    ## Y and X are such that we regress to get gamma
    Y <- -as.vector(ft.mat)
    X <- as.vector(fa.mat)
##     W <- as.vector(fat)

    if (is.null(smoothed))
    {

        m <- lm(Y ~ -1 + coh:X)
##        m <- lm(Y ~ -1 + coh:X, weights = W)
##         print("est weighted")
        coef.mat <- summary(m)$coef
        my.coef <- coef.mat[,"Estimate"]
        my.se <- coef.mat[,"Std. Error"]
        mc <- gsub(":X", "", names(my.coef))
        mc <- gsub("coh", "", mc)

        gamma <- my.coef
    }
    if (!is.null(smoothed))
    {
        se.smooth <- gamma.smooth <- rep(NA, length(coh.names))
        for (i in seq(coh.names))
        {
            if(i > 2 && i < length(coh.names)-2)
            {
                pre.last.cohort <- coh.names[i-2]
                last.cohort <- coh.names[i-1]
                this.cohort <- coh.names[i]
                next.cohort <- coh.names[i+1]
                post.next.cohort <- coh.names[i+2]


                wt <- rep(0, length(X))
                ## wt[coh == last.cohort] <- .2
                ## wt[coh == pre.last.cohort] <- .2
                ## wt[coh == this.cohort] <- .2
                ## wt[coh == next.cohort] <- .2
                ## wt[coh == post.next.cohort] <- .2
                wt[coh == last.cohort] <- 1/3
                wt[coh == this.cohort] <- 1/3
                wt[coh == next.cohort] <- 1/3

                m <- lm(Y ~ -1 + X, weights = wt)


                gamma.smooth[i] <- coef(m)
                if (!is.na(coef(m)))
                    se.smooth[i] <- summary(m)$coef[,"Std. Error"]
            }
        }
        gamma <- gamma.smooth
        my.se <- se.smooth
    }

    ## cap gamma at 1/3, which caps Sp at 1/2
    if (any(gamma[!is.na(gamma)] > 1/3))
        print("Warning: very big values of Sp, > 1/2")
    ##     gamma[gamma > 1/3] <- 1/3
    ##     gamma[gamma < -1/3] <- -1/3

    ## convert back to Sp
    Sp <- gamma/(1-gamma)

    ## make sure we return full length
    if (is.null(smoothed))
    {
        Sp.vec <- rep(NA, length(coh.names))
        Sp.vec[coh.names %in% mc] <- Sp
        names(Sp.vec) <- coh.names
        Sp.se.vec <- rep(NA, length(coh.names))
        Sp.se.vec[coh.names %in% mc] <- my.se
        names(Sp.se.vec) <- coh.names
        ## get rid of NAs, replacing with zero
        Sp.vec[is.na(Sp.vec)] <- 0
    }
    if (!is.null(smoothed))
    {
        Sp.vec <- Sp
        names(Sp.vec) <- coh.names
        Sp.se.vec <- my.se
        names(Sp.se.vec) <- coh.names
    }

    if (return.se)
        return(Sp.se.vec)
    if (!return.se)
        return(Sp.vec)
}

get.Sp.vec.truncated.means <- function(fat, eps = .01) {
    ## output: Sp.vec should have length nrow(Fat) + ncol(Fat) - 1,
    ## just as if we went through diagonals ...

    ## fat =  matrix of fertility rates by age a and time t, with
    ## numeric labels.

    ## eps = tolerance for numerator and denominator of estimation
    ## formula.  if fx, next.fx, or last.fx too small assign Sp <- 0


    ## get cohort matrix set up, as well as other useful vectors
    fac.big <- per2coh(fat)
    fac <- fac.big[,!apply(fac.big, 2, all.na)] # shrink to right size
    cohort <- colnames(fac)
    x <- 1:nrow(fac)                    # age, arbitrary
    Sp.vec <- rep(NA, length(cohort))
    names(Sp.vec) <- cohort

    ## now loop through cohort matrix
    for (i in 1:length(cohort))
    {
        if (i != 1 && i != length(cohort)) {
            ##       i <- 20
            fx <- fac[,i]
            next.fx <- fac[,i+1]
            last.fx <- fac[,i-1]


            h <- max(x[!is.na(next.fx)])
            l <- min(x[!is.na(last.fx)])
##            h <- max(x[!is.na(fx)])
##            l <- min(x[!is.na(fx)])

            s <- x %in% l:h

            mu <- sum((x * fx)[s])/ sum(fx[s])
            next.mu <- sum((x * next.fx)[s])/ sum(next.fx[s])
            last.mu <- sum((x * last.fx)[s])/ sum(last.fx[s])
            mu.prime <- (next.mu - last.mu)/2

            vhc <- fx[x %in% h]/sum(fx[s])
            vlc <- fx[x %in% l]/sum(fx[s])

            adjust <- 1 + vhc*(mu-h) + vlc*(l-mu)
            Sp <- mu.prime / adjust

            ## if fx, next.fx, or last.fx too small
            ## assign Sp <- 0
            if (sum(fx, na.rm = T) < eps ||
                sum(next.fx, na.rm = T) < eps ||
                sum(last.fx, na.rm = T) < eps)
                Sp <- 0
            ## save current value
            Sp.vec[i] <- Sp

        }
    }
    ## replace any NaN or NA with zero
    Sp.vec[is.na(Sp.vec)] <- 0
    Sp.vec[is.nan(Sp.vec)] <- 0
    return(Sp.vec)
}

get.Sp.vec <- get.Sp.vec.truncated.means


vec2diagmat <- function(vec, mat) {
    ## we want to take vec 1:3
    ## and turn it into
    ## 2 3
    ## 1 2
    if (length(vec) != nrow(mat) + ncol(mat) -1) {
        print("error: vec and mat not conformable")
        return(NULL);
    }
    D.gen <- col(mat) - row(mat)
    D <- D.gen - min(D.gen) + 1
    diagmat <- matrix(NA, nrow(mat), ncol(mat))
    for (i in 1:nrow(mat)) {
        diagmat[i,] <- vec[D[i,]]
    }
    diagmat
}




get.tfr.dagger <- function(fat, n.iter = 5,
                           get.Sp.fun = get.Sp.vec.truncated.means, ...) {
    tfr.dagger.mat <- matrix(NA, n.iter, ncol(fat))
    Sp.vec.mat <- matrix(NA, n.iter, ncol(fat) + nrow(fat) -1)
    Sp.se.mat <- matrix(NA, n.iter, ncol(fat) + nrow(fat) -1)
    Sp.diag.mat <- matrix(NA, nrow(fat), ncol(fat))

    ## get Sp.vec
##    Sp.vec <- get.Sp.vec(fat)
    ## put in a diagonal matrix in order to then calculate tfr dagger
##    Sp.diag.mat <- vec2diagmat(Sp.vec, fat)



    ## Iterate qt and S.prime
    tfr <- apply(fat, 2, sum)
    tfr.dag  <- tfr                     # initialize
    fat.tilde <- fat                    # initialize
    for (iter in 1:n.iter)
    {
        ## assign results of last iteration
        Sp.vec <- get.Sp.fun(fat.tilde, ...)
        Sp.vec.mat[iter,] <- Sp.vec
##        Sp.se.mat[iter,] <- Sp.se
        Sp.diag.mat <- vec2diagmat(Sp.vec, fat)
        for (i in 1:ncol(fat))
        {
            tfr.dag[i] <- sum(fat[,i] * (1 + Sp.diag.mat[,i]))
        }
        tfr.dagger.mat[iter,] <- tfr.dag


        qt.est <- tfr.dag
        ## now fill in missing values by extending
        qt.est <- approx(seq(qt.est), qt.est, seq(qt.est), rule = 2)$y

        for (i in 1:ncol(fat)) {
            fat.tilde[,i] <- fat[,i]/qt.est[i]
        }

        ## where qt is na, don't change fat.tilde
##         fat.tilde[is.na(fat.tilde)] <- fat[is.na(fat.tilde)]
##        print("dim(fat.tilde)")
##        print(dim(fat.tilde))

##         Sp.vec <- get.Sp.vec(fat.tilde, ...)
        ## defaults to get.Sp.vec.truncated.means()
##         Sp.vec <- get.Sp.fun(fat.tilde, ...)

##        print(paste("length(Sp.vec)", length(Sp.vec)))
##        print(head(Sp.vec))
##        print(tail(Sp.vec))

        ## assign se if we are using get.Sp.vec with regression
##         ifelse(!is.null(formals(get.Sp.vec, ...)$return.se),
##               Sp.se <- NULL,
        ##               Sp.se <- get.Sp.vec(fat.tilde, return.se = TRUE))
        ## Sp.se <- NULL
        ## Sp.se <- rep(NA, ncol(Sp.se.mat))
        ## Sp.vec.mat[iter,] <- Sp.vec
        ## Sp.se.mat[iter,] <- Sp.se
        ## Sp.diag.mat <- vec2diagmat(Sp.vec, fat)

        ## for (i in 1:ncol(fat))
        ## {
        ##     tfr.dag[i] <- sum(fat[,i] * (1 + Sp.diag.mat[,i]))
        ## }
        ## tfr.dagger.mat[iter,] <- tfr.dag
    }
    return(list(tfr.dagger = tfr.dagger.mat[n.iter,],
                tfr.dagger.mat = tfr.dagger.mat,
                Sp.vec.mat = Sp.vec.mat,
                Sp.se.mat = Sp.se.mat,
                Sp.vec = Sp.vec))
}


get.tfr.dagger.orig <- function(fat, n.iter = 5) {
    ##
##    fat <- fat
##    n.iter <- 5
##    print("dim(fat)")
##    print(dim(fat))
    tfr.dagger.mat <- matrix(NA, n.iter, ncol(fat))
    Sp.vec.mat <- matrix(NA, n.iter, ncol(fat) + nrow(fat) -1)
    Sp.se.mat <- matrix(NA, n.iter, ncol(fat) + nrow(fat) -1)
    Sp.diag.mat <- matrix(NA, nrow(fat), ncol(fat))
    ## get Sp.vec
##    Sp.vec <- get.Sp.vec(fat)
    ## put in a diagonal matrix in order to then calculate tfr dagger
##    Sp.diag.mat <- vec2diagmat(Sp.vec, fat)



    tfr <- apply(fat, 2, sum)
    tfr.dag  <- tfr                     # initialize
    fat.tilde <- fat                    # initialize
    for (iter in 1:n.iter)
    {
        qt.est <- tfr.dag
        ## now fill in missing values by extending

##         qt.est <- approx(seq(qt.est), qt.est, seq(qt.est), rule = 2)$y

        for (i in 1:ncol(fat))
            fat.tilde[,i] <- fat[,i]/qt.est[i]
##        print("dim(fat.tilde)")
##        print(dim(fat.tilde))

        Sp.vec <- get.Sp.vec.orig(fat.tilde)

        Sp.vec.mat[iter,] <- Sp.vec
        Sp.diag.mat <- vec2diagmat(Sp.vec, fat)

        for (i in 1:ncol(fat))
        {
            tfr.dag[i] <- sum(fat[,i] * (1 + Sp.diag.mat[,i]))
        }
        tfr.dagger.mat[iter,] <- tfr.dag
    }
    return(list(tfr.dagger = tfr.dagger.mat[n.iter,],
                tfr.dagger.mat = tfr.dagger.mat,
                Sp.vec.mat = Sp.vec.mat,
                Sp.vec = Sp.vec))
}





shift.mat <- function(x, fx.mat, Kt, fat, fac)
{
    ## x is age
    ## fx.mat is a matrix of identical f0 for every year
    ## Kt is a generalized notation for shifts (the Sc, the Rt or the Sc and Rt, if you will)
    ## fat and fac are matrices to provide the right dimnames in the
    ## output

    ## convert fx.mat to long vector form
    fx.vec <- as.vector(fx.mat)
    ## replace zeros with something very small so that we can take logs
    eps <- 10^-10
    fx.vec[fx.vec == 0] <- eps

    ## convert x to long vector form
    x.vec <- x[row(fx.mat)]

    ## now get Kt in long vector for
    ## check errors in length of Kt, which can either be as long as
    ## the fx.mat or as long as there are columns in fx.mat (e.g.,
    ## with BF model)
    stopifnot(length(Kt) == length(fx.mat) ||
              length(Kt) == ncol(fx.mat))
    if (length(Kt) == length(fx.mat))   # if Kt is matrix
        Kt.vec <- as.vector(Kt)
    if (length(Kt) == ncol(fx.mat)) {
        ## next is probably right but tricky
        ## we want
        ## Kt[1] Kt[2] Kt[3]  ...
        ## Kt[1] Kt[2] Kt[3]  ...
        ## ...
        Kt.vec <- Kt[col(fx.mat)]
    }
    ## Originally we had rule = 2, but Tom was worried
    ##     fxout.vec <- approx(x.vec, fx.vec, xout = x.vec - Kt.vec,
    ##                        rule = 2)$y

    ## Now we have splines:
    ## Here we follow r help page where it says natural method is good
    ## for extrapolation -- that it uses linear extrapolation, which
    ## we make log-linear by working in logs.  note further that the
    ## spline does not act like oldest age and youngest age are next
    ## to each other, since x.vec goes from oldest age to youngest
    ## age, and spline knows that these are distant by e.g. 54 to
    ## 14. In other words it uses the actual ages, not just the
    ## indices of the vector. BUT DOES IT USE ADJACENT AGES FOR
    ## COHORTS AND SUCH?

    ##    print(log(fx.vec))
    if (sum(is.infinite(fx.vec)) > 0) {
        print("error in fx.vec: infinite values")
        print((fx.vec))
    }

    ## for debugging
    verbose = FALSE
    if (verbose) {
        print(length(x.vec))
        print(Kt.vec)
        print(cbind(x.vec, log(fx.vec), x.vec - Kt.vec))
      }
    fxout.vec <- exp(spline(x = x.vec, y = log(fx.vec), n = length(x.vec),
                       xout = x.vec - Kt.vec, method = "natural")$y)

    if (0) { ## run this to verify that splining with ages is continuous
    x = c(0:5, 0:5)
    xx <- c(seq(0,5, .1),seq(0,5, .1))
    t.out = spline(x = x, y = x^3,
    xout = xx - c(rep(.1, length(xx)/2),rep(.1, length(xx)/2)),
    method = "natural")
    y.out = t.out$y
    plot(x, x^3, type = "p")
    lines(t.out, type = "o", cex = .5, pch = 19, col = "red")
}

    ## old:
    ## here's a rule = 1, (NA) alternative, that is modified to give zeros
    ##    fxout.vec <- approx(x.vec, fx.vec, xout = x.vec - Kt.vec,,
    ##                        rule = 1)$y
    ##  fxout.vec[is.na(fxout.vec)] <- 0

    fxout.mat <- matrix(fxout.vec, nrow(fx.mat), ncol(fx.mat))
    if (all(dim(fxout.mat) == dim(fat)))
        dimnames(fxout.mat) <- dimnames(fat)
    if (all(dim(fxout.mat) == dim(fac)))
        dimnames(fxout.mat) <- dimnames(fac)

    fxout.mat
}

opt.fun.2.free <- function(p, f0, Sc, Rt, fat, fac)
{

    ## idea here is to set one parameter at at time free (in addition
    ## of course to qt). we figure out which parameter is free, and
    ## then assign p appropriately. then we get fat.hat, and then we
    ## return the objective


    ## usage: opt.fun.1.free(p, f0 = NULL, Sc, Rt)
    if(is.null(f0)) {
        f0 = antilogit(p[1:nrow(fat)])
        qt = p[(nrow(fat)+1):length(p)]
    }
    if(is.null(Sc)) {
        Sc = p[1:ncol(fac)]
        qt = p[(ncol(fac)+1):length(p)]
    }
    if(is.null(Rt)) {
        Rt = p[1:ncol(fat)]
        qt = p[(ncol(fat)+1):length(p)]
    }

    fat.hat <- p2fat.hat(fat=fat, fac=fac,
                         f0 = f0,
                         Rt = Rt,
                         Sc = Sc,
                         qt = qt)
    obj <- sum((fat-fat.hat)^2, na.rm = T)
    obj
}

opt.fun.1.free <- function(p, f0, Sc, Rt, qt, fat, fac)
{

    ## qt can also be free

    # we figure out which parameter is free, and
    ## then assign p appropriately. then we get fat.hat, and then we
    ## return the objective


    ## usage: opt.fun.1.free(p, f0 = NULL, Sc, Rt)
    if(is.null(f0)) {
        f0 = antilogit(p)
    }
    if(is.null(Sc)) {
        Sc = p
    }
    if(is.null(Rt)) {
        Rt = p
    }
    if(is.null(qt)) {
        qt = p
    }

    fat.hat <- p2fat.hat(fat=fat, fac=fac,
                         f0 = f0,
                         Rt = Rt,
                         Sc = Sc,
                         qt = qt)
##     pen <- sum(diff(Sc)^2, na.rm = T) * .3
##        pen <- sum(diff(diff(Sc))^2, na.rm = T) * .03
    pen <- 0

    obj <- sum((fat-fat.hat)^2, na.rm = T) + pen
    obj
}

opt.fun.norm.1.free <-
function(p, f0.param, Sc, Rt, qt, fat, fac)
{
    # figure out which parameter is free, then assign p
    # appropriately. then get fat.hat, and then return objective

    ## usage: opt.fun.1.free(p, f0 = NULL, Sc, Rt)
    if(is.null(f0.param)) {
        mu <- p[1]
        sigma <- p[2]
        x <- as.numeric(rownames(fat))
        f0 <- dnorm(x, mean = mu, sd = sigma)
        f0.param <- p
##        f0 = antilogit(p)
    }
    if(is.null(Sc)) {
        Sc = p
    }
    if(is.null(Rt)) {
        Rt = p
    }
    if(is.null(qt)) {
        qt = p
    }


    mu <- f0.param[1]
    sigma <- f0.param[2]

    f0 <- dnorm(x, mean = mu, sd = sigma)

    fat.hat <- p2fat.hat(fat=fat, fac=fac,
                         f0 = f0,
                         Rt = Rt,
                         Sc = Sc,
                         qt = qt)
##     pen <- sum(diff(Sc)^2, na.rm = T) * .3
##        pen <- sum(diff(diff(Sc))^2, na.rm = T) * .03
    pen <- 0

    obj <- sum((fat-fat.hat)^2, na.rm = T) + pen
    obj
}

opt.fun.1.free.gx <- function(p, f0, Sc, Rt, qt, gx, fat, fac)
{

    ## qt can also be free

    # we figure out which parameter is free, and
    ## then assign p appropriately. then we get fat.hat, and then we
    ## return the objective


    ## usage: opt.fun.1.free(p, f0 = NULL, Sc, Rt)
    if(is.null(f0)) {
        f0 = antilogit(p)
    }
    if(is.null(Sc)) {
        Sc = p
    }
    if(is.null(Rt)) {
        Rt = p
    }
    if(is.null(qt)) {
        qt = p
    }
    if(is.null(gx)) {
        gx = p
        ## freeze gx at 1 for ages 24 to 38
##         age = rownames(fat)
##        s = age %in% 24:38
##        gx[s] <- 1
    }


    fat.hat <- p2fat.hat.gx(fat=fat, fac=fac,
                            f0 = f0,
                            Rt = Rt,
                            Sc = Sc,
                            qt = qt,
                            gx = gx)
    obj <- sum((fat-fat.hat)^2, na.rm = T)
    obj
}


per2c <- function(A)
{
    col(A) - row(A) + nrow(A)
}

p2fat.hat <- function(fat, fac, f0, Rt, Sc, qt)
{
    x <- 1:length(f0)
    rt <- diff(c(0,Rt))
    Rat.vec <- Rt[col(fat)]
##    print(length(Rat.vec))
    Sat.vec <- Sc[per2c(fat)]           # note per2c returns diag number
##    print(length(Sat.vec))
    Rat.plus.Sat <- Rat.vec + Sat.vec
##    print(f0)
##    print(row(fat))
    if (sum(is.infinite(f0)) > 0) {
        print("problem with infinite values in f0")
        print((f0))
    }

    fx.mat <- matrix(f0[row(fat)], nrow(fat), ncol(fat))
##    if (sum(is.infinite(fx.mat)) > 0)
##        print((fx.mat))

##    print(dim(fx.mat))
    fat.shift <- shift.mat(x, fx.mat, Kt=Rat.plus.Sat, fac=fac, fat=fat)
    qat <- qt[col(fat)]
    rat <- rt[col(fat)]

    ##     fat.hat <- t(t(fat.shift) * qt * (1-rt)) # t(t()) makes multiplicationra
    fat.hat <- fat.shift * qat * (1-rat)
    fat.hat
}

p2fat.hat.gx <- function(fat, fac, f0, Rt, Sc, qt, gx)
{
    x <- 1:length(f0)
    rt <- diff(c(0,Rt))
    Rat.vec <- Rt[col(fat)]
##    print(length(Rat.vec))
    Sat.vec <- Sc[per2c(fat)]
##    print(length(Sat.vec))
    Rat.plus.Sat <- Rat.vec + Sat.vec
##    print(f0)
##    print(row(fat))
    fx.mat <- matrix(f0[row(fat)], nrow(fat), ncol(fat))
##    print(dim(fx.mat))
    fat.shift <- shift.mat(x, fx.mat, Rat.plus.Sat, fac=fac, fat=fat)
    qat <- qt[col(fat)]
    rat <- rt[col(fat)]
    gxt <- gx[row(fat)]

    ##     fat.hat <- t(t(fat.shift) * qt * (1-rt)) # t(t()) makes multiplicationra
    fat.hat <- fat.shift * qat * (1-rat) *gxt
    fat.hat
}

opt.bf <- function(fat, ...)
{
    print("doing opt.bf")
    ## get fac
    fac.big <- per2coh(fat)
    fac <- fac.big[,!apply(fac.big, 2, all.na)] # shrink to right size
    print("bango")
    ## starting values
    obj.last = Inf
    ## fix many
    bf.out <- bf.fit.simple(fat)
    rt.hat <- bf.out$rt
    ## replace any NA rt's with zero
    rt.hat[is.na(rt.hat)] <- 0
    Rt.hat <- cumsum(rt.hat)                # bf fit
    qt.hat <- bf.out$tfr.star
    ## set Sc.hat at zero, since fitting BF
    Sc.hat <- rep(0, ncol(fac))
    ## initialize f0, which will be free
    f0.hat <- apply(fat, 1, mean)           # mean by age
  eps <- 10^-10
    f0.hat <- f0.hat + eps
    ## f0
    out <- nlm(f0 = NULL, Sc = Sc.hat, Rt = Rt.hat,
               qt = qt.hat,
               f = opt.fun.1.free, fac = fac, fat = fat,
               p = logit(f0.hat), ...)
    f0.hat <- antilogit(out$est)
##    print(paste(round(out$m, 5), i, "f0 step"))
    fat.hat <- p2fat.hat(fat, fac, f0 = f0.hat,
                         Rt = Rt.hat,
                         Sc = Sc.hat,
                         qt = qt.hat)

    ## return results
    result <- list("out" = out,
                   "obj" = out$m,
                   "f0" = f0.hat,
                   "fat.hat" = fat.hat)
}

opt.gc <- function(fat, ...)
{
    print("doing opt.gc")
    ## get fac
    fac.big <- per2coh(fat)
    fac <- fac.big[,!apply(fac.big, 2, all.na)] # shrink to right size

    ## starting values
    ## fix many
    gc.out <- get.tfr.dagger(fat)
    Sp.vec <- get.tfr.dagger(fat)$Sp.vec    # gc fit
    Sp.vec[is.na(Sp.vec)] <- 0
    ## replace any NaN Sp's with zero
    Sp.vec[is.nan(Sp.vec)] <- 0
    Sp.vec[is.infinite(Sp.vec)] <- 0


    Sc.hat <- cumsum(Sp.vec)
    qt.hat <- gc.out$tfr.dagger

    ## set Rt.hat at zero, since fitting GC
    Rt.hat <- rep(0, ncol(fat))
    ## initialize f0, which will be free
    f0.hat <- apply(fat, 1, mean)           # mean by age
    eps <- 10^-6
    f0.hat <- f0.hat + eps
    ## f0

    ## not the problem
    ## if (!complete.cases(Sc.hat))
##     {
##         print("!complete.cases(Sc.hat))")
##        print(Sc.hat)
##    }
##    Sc.hat.tmp <<- Sc.hat

    out <- nlm(f0 = NULL, Sc = Sc.hat, Rt = Rt.hat,
               qt = qt.hat,
               f = opt.fun.1.free, fac = fac, fat = fat,
               p = logit(f0.hat), ...)
    f0.hat <- antilogit(out$est)
##     print(paste(round(out$m, 5), i, "f0 step"))
    fat.hat <- p2fat.hat(fat, fac, f0 = f0.hat,
                         Rt = Rt.hat,
                         Sc = Sc.hat,
                         qt = qt.hat)

    ## return results
    result <- list("out" = out,
                   "obj" = out$m,
                   "f0" = f0.hat,
                   "fat.hat" = fat.hat)

}

opt.gc.full <- function(fat, maxiter = 10, mpl = 1, ...)
{
    print("doing opt.gc.full")
    ## get fac
    fac.big <- per2coh(fat)
    fac <- fac.big[,!apply(fac.big, 2, all.na)] # shrink to right size

    ## starting values
    ## fix many
    gc.out <- get.tfr.dagger(fat)
    Sp.vec <- get.tfr.dagger(fat)$Sp.vec    # gc fit
    Sp.vec[is.na(Sp.vec)] <- 0
    ## replace any NaN Sp's with zero
    Sp.vec[is.nan(Sp.vec)] <- 0
    Sp.vec[is.infinite(Sp.vec)] <- 0


    Sc.hat <- cumsum(Sp.vec)
    qt.hat <- gc.out$tfr.dagger

    ## set Rt.hat at zero, since fitting GC
    Rt.hat <- rep(0, ncol(fat))
    ## initialize f0, which will be free
    f0.hat <- apply(fat, 1, mean)           # mean by age
    eps <- 10^-6
    f0.hat <- f0.hat + eps
    ## f0

    ## not the problem
    ## if (!complete.cases(Sc.hat))
##     {
##         print("!complete.cases(Sc.hat))")
##        print(Sc.hat)
##    }
##    Sc.hat.tmp <<- Sc.hat


    for (j in 1:maxiter) {

        ## f0
        out <- nlm(f0 = NULL, Sc = Sc.hat, Rt = Rt.hat,
                   qt = qt.hat,
                   f = opt.fun.1.free, fac = fac, fat = fat,
                   p = logit(f0.hat), print.level = mpl)
        f0.hat <- antilogit(out$est)
        print(paste(round(out$m, 5), j, "f0 step"))

        ## Sc
        out <- nlm(f0 = f0.hat, Sc = NULL, Rt = Rt.hat,
                   qt = qt.hat,
                   f = opt.fun.1.free, fac = fac, fat = fat,
                   p = Sc.hat, print.level = mpl)
        Sc.hat <- out$est
        print(paste(round(out$m, 5), j, "Sc step"))


        ## f0
        out <- nlm(f0 = NULL, Sc = Sc.hat, Rt = Rt.hat,
                   qt = qt.hat,
                   f = opt.fun.1.free, fac = fac, fat = fat,
                   p = logit(f0.hat), print.level = mpl)
        f0.hat <- antilogit(out$est)
        print(paste(round(out$m, 5), j, "f0 step"))

        ## ## Rt
        ## out <- nlm(f0 = f0.hat, Sc = Sc.hat, Rt = NULL,
        ##            qt = qt.hat,
        ##            f = opt.fun.1.free, fac = fac, fat = fat,
        ##            p = Rt.hat, print.level = mpl)
        ## Rt.hat <- out$est
        ## print(paste(round(out$m, 5), j, "Rt step"))

        ## qt
        out <- nlm(f0 = f0.hat, Sc = Sc.hat, Rt = Rt.hat,
                   qt = NULL,
                   f = opt.fun.1.free, fac = fac, fat = fat,
                   p = qt.hat, print.level = mpl)
        qt.hat <- out$est
        print(paste(round(out$m, 5), j, "qt step"))

    }

    ## out <- nlm(f0 = NULL, Sc = Sc.hat, Rt = Rt.hat,
    ##            qt = qt.hat,
    ##            f = opt.fun.1.free, fac = fac, fat = fat,
    ##            p = logit(f0.hat), ...)
    ## f0.hat <- antilogit(out$est)
##     print(paste(round(out$m, 5), i, "f0 step"))
    fat.hat <- p2fat.hat(fat, fac,
                         f0 = f0.hat,
                         Rt = Rt.hat,
                         Sc = Sc.hat,
                         qt = qt.hat)

    ## return results
    result <- list("out" = out,
                   "obj" = out$m,
                   "f0" = f0.hat,
                   "Rt" = Rt.hat,
                   "Sc" = Sc.hat,
                   "qt" = qt.hat,
                   "fat.hat" = fat.hat)

}
