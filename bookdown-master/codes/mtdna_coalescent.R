## ----------------------------------------------------------------------------------------------------------------------------
source("read_mtdna.R")
haps <- seqs


## ----------------------------------------------------------------------------------------------------------------------------
print(nchar(haps[213]))
a_segment = substr(haps[213], 1, 100)
print(a_segment)
## tip: don't try to print the 16,000 character whole string. it will clog up your computer.


## ----------------------------------------------------------------------------------------------------------------------------
my.list <- strsplit(haps, "")
H <- do.call(cbind, my.list)
print(H[1:10, 1:4]) ## all the same


## ----------------------------------------------------------------------------------------------------------------------------
hap1 = H[,1]
hap2 = H[,2]
s <- min(which(hap1 != hap2))
head(H[s + -2:2, 1:4]) ## the polymorphic site in context


## ----------------------------------------------------------------------------------------------------------------------------
table(H[s,])


## ----------------------------------------------------------------------------------------------------------------------------

set.seed(1)
hap_ids = 1:ncol(H)
hap_id_sample = sample(hap_ids,
                       size = 200,
                       replace = FALSE)

hap_id.mat <- matrix(hap_id_sample, 100, 2)

pairwise_diff_fun <- function(hap1, hap2)
{
    h1 <- hap1
    h2 <- hap2
##    h1 <- unlist(strsplit(hap1, ""))
##    h2 <- unlist(strsplit(hap2, ""))
    h1[h1 == "N"] <- NA ## note "N" means missing
    h2[h2 == "N"] <- NA ## making these NA avoids counting as polymorphism
    k = sum(h1 != h2, na.rm = T)
    n_valid = sum(!is.na(h1) & !is.na(h2))
    return(list(k = k, n_valid = n_valid))
}

pairwise_diff_fun(H[,1], H[,2])


## ----------------------------------------------------------------------------------------------------------------------------

P.vec = NULL
C.vec = NULL
for (i in 1:nrow(hap_id.mat))
{
    hap_id.1 = hap_id.mat[i,1]
    hap_id.2 = hap_id.mat[i,2]
    hap1 = H[,hap_id.1]
    hap2 = H[,hap_id.2]
    out = pairwise_diff_fun(hap1, hap2)
    P.vec[i] = out$k
    C.vec[i] = out$n_valid
}

Y.bar = P.vec/C.vec
head(P.vec)
head(C.vec)
head(Y.bar)


## ----------------------------------------------------------------------------------------------------------------------------
theta_m = 2.21 * 10^(-8) ## Batini page 6 (TMRCA estimation)
T.vec <- -(1/2) * (1/theta_m) * log(1 - Y.bar) ## TMRCAs in years ago
head(T.vec, n = 10)


## ----------------------------------------------------------------------------------------------------------------------------
hist(T.vec)


## ----------------------------------------------------------------------------------------------------------------------------
## Plot survival curve by order of T
St = (100:1)/100 ## or more generally (length(T.vec):1)/length(T.vec)
t = kya = sort(T.vec)/1000
plot(kya, St, type = "l",
     xlab = "Kilo years ago", ylab = "Fraction of pairs without common ancestor",
     main = "Estimated probability of not coalescing")
plot(kya, log(St), type = "l",
     xlab = "Kilo years ago", ylab = "Log fraction of pairs without common ancestor",
     main = "Etimatated probability of not coalescing, log scale")
     


## ----------------------------------------------------------------------------------------------------------------------------
out = lowess(x = kya, y = St, f = 1/5)

St_smooth = out$y
kya_smooth = out$x
plot(kya, St, cex = .5)
lines(kya_smooth, St_smooth, type = 'l')


## ----------------------------------------------------------------------------------------------------------------------------
haz_hat = -diff(St_smooth)/diff(kya_smooth)
plot(kya_smooth[-1], haz_hat, type = 'l')


## ----------------------------------------------------------------------------------------------------------------------------
## we choose these time boundaries arbitrarily ... not sure if
## we'll be able to see the "expansion" after ice age ...
## x = c(0, 12,  20, 40, 65,  180 )*1000    # time interval boundaries
x = c(0,2, 5, 10,  20, 30, 40, 65, 180) * 1000 ## time interval boundaries


## ----------------------------------------------------------------------------------------------------------------------------
get_nax <- function(Ti, x)
{
    ## get person years lived in interval by those who die
    nax <- NULL
    for (i in 1:(length(x)-1))
    {
        s <- Ti >= x[i] & Ti < x[i+1]
        nax[i] = mean(Ti[s] - x[i])
    }
    return(nax)
}


## ----------------------------------------------------------------------------------------------------------------------------
n <- diff(x)
T.vec.by.cat <- cut(T.vec, x, include.lowest = T, right = F)
ndx = table(T.vec.by.cat)
lx = rev(cumsum(rev(ndx)))
lxpn = c(lx[-1], 0)
nax = get_nax(Ti = T.vec, x = x)
nLx = n*lxpn + nax * ndx ## exposure
nmx = ndx/nLx ## hazard
lt <- cbind(x = x[-length(x)], n, ndx, lx, nax, nLx, nmx)
print(lt)


## ----------------------------------------------------------------------------------------------------------------------------
x.mid = x[-length(x)] + n/2
plot(x.mid, nmx, type = 'o')
axis(2)
lines(kya_smooth[-1] * 1000, haz_hat/1000, type = "l")


## ----------------------------------------------------------------------------------------------------------------------------

Ne_smooth = 2 / (haz_hat/1000 * 25)

Ne_lifetable = 2 / (nmx * 25)

## create step function for plotting
Ne_lifetable_step = rep(Ne_lifetable, n/1000)
kya_step = 1:(max(x)/1000)



## ----------------------------------------------------------------------------------------------------------------------------
plot(kya_smooth[-1], Ne_smooth, type = 'l', ylim = c(1000, 60000), log = 'y',
     lty = 2)
lines(kya_step, Ne_lifetable_step, type = 'l', lwd = 2)


## ----------------------------------------------------------------------------------------------------------------------------
n_trials = 40
Ne.mat <- matrix(NA, nrow = n_trials, ncol = max(x)/1000)
set.seed(1)
for (r in 1:n_trials)
{
    print(r)
    ## sample
    hap_id_sample = sample(hap_ids,
                           size = 200,
                           replace = FALSE)
    hap_id.mat <- matrix(hap_id_sample, 100, 2)
    ## estimate MRCA distribution
    P.vec = NULL
    C.vec = NULL
    for (i in 1:nrow(hap_id.mat))
    {
        hap_id.1 = hap_id.mat[i,1]
        hap_id.2 = hap_id.mat[i,2]
        hap1 = H[,hap_id.1]
        hap2 = H[,hap_id.2]
        out = pairwise_diff_fun(hap1, hap2)
        P.vec[i] = out$k
        C.vec[i] = out$n_valid
    }
    Y.bar = P.vec/C.vec
    T.vec <- -(1/2) * (1/theta_m) * log(1 - Y.bar) ## TMRCAs in years ago
    ## estimate Ne
    T.vec.by.cat <- cut(T.vec, x, include.lowest = T, right = F)
    ndx = table(T.vec.by.cat)
    lx = rev(cumsum(rev(ndx)))
    lxpn = c(lx[-1], 0)
    nax = get_nax(Ti = T.vec, x = x)
    nLx = n*lxpn + nax * ndx ## exposure
    nmx = ndx/nLx ## hazard
    Ne_lifetable = 2 / (nmx * 25)
    ## create step function for plotting
    Ne_lifetable_step = rep(Ne_lifetable, n/1000)
    kya_step = 1:(max(x)/1000)
    ## save result
    Ne.mat[r,] <- Ne_lifetable_step
}

Ne.interval <- apply(Ne.mat, 2, quantile, c(.1,.5, .9))
matplot(t(Ne.interval), type = 'l', log = 'y', col = "grey", lty = 1, lwd = 2)
lines(Ne.interval["50%",], lwd = 4)

