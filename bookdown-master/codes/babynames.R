## Fisher-Wright simulation of Baby Names (Hahn and Bentley)

###############
## data prep ##
###############


download.file(url= "https://www.ssa.gov/oact/babynames/names.zip",
              "./names.zip")
unzip("names.zip", exdir = "./names")

library(data.table)

filenames <- system("ls ./names/*.txt", intern = T)
mylist <- vector("list", length(filenames))
names(mylist) <- gsub(pattern = "[^0-9]", replace = "", filenames)
for (i in 1:length(filenames))
{
    myfile <- filenames[i]
    mylist[[i]] <- fread(myfile)
}

dt <- rbindlist(mylist, idcol = "year")
names(dt) <- c("year", "name", "sex", "N")

## dt
##          year      name sex    N
##       1: 1880      Mary   F 7065
##       2: 1880      Anna   F 2604
##       3: 1880      Emma   F 2003
##       4: 1880 Elizabeth   F 1939
##       5: 1880    Minnie   F 1746
##      ---
## 1957042: 2018     Zylas   M    5
## 1957043: 2018     Zyran   M    5
## 1957044: 2018     Zyrie   M    5
## 1957045: 2018     Zyron   M    5
## 1957046: 2018     Zzyzx   M    5

## ok, we have the data now

###############################
## Plot observed frequencies ##
###############################

## male 1900-1909

my.dt <- dt[sex == "M" & year %in% 1900:1909]
foo <- my.dt[, .(N = sum(N)), by = name]
foo <- foo[order(N, decreasing = T)]
bar <- foo[1:1000,] ## 1000 top names

## now let's do a power law plot

my.breaks <- c(0, 2^(0:11)/10000)
bar[, p := round(prop.table(N),5)]
bar[, pcat := cut(p, breaks =  my.breaks, right = F, include.lowest = T)]

out <- unclass(prop.table(table(bar$pcat)))


my.x <- my.breaks[-length(my.breaks)] + diff(my.breaks)/2
plot(my.x, out, log = "xy")



###################
## FW simulation ##
###################


mut <- function(x, mu)
{
    ## m, the individuals that mutate
    m <- which(rbinom(length(x), 1, mu) == 1)
    if (length(m) == 0)
        return(x)
    suffix <- 10000*round(runif(length(m)),4)
    x[m] <- paste0(x[m], ".", suffix)
    x
}

fwm <- function(N, n_gen, mu = 0)
{
x <- paste(1:N)
A <- matrix(NA, nrow = n_gen, ncol = N)
for (i in 1:n_gen)
{
    A[i,] <- x
    x <- sample(x, size = N, replace = T)
    x <- mut(x, mu)
    x
}
return(A)
}

## let's look at evolution over time of G:
## chance that two individuals are of same type

get.G <- function(x)
{
    tt <- table(x)
    p <- prop.table(tt)
    sum(p^2)
}


## without mutation
A <- fwm(1000, n_gen = 4000, mu = 0)
G.vec <- apply(A, 1, get.G)
plot(G.vec)

## with mutation, 1 trial

N = 1000
A <- fwm(N, n_gen = 3000, mu = 4/N)
G.vec <- apply(A, 1, get.G)
plot(G.vec)

## average over 100 trials

n_gen = 2000
n_trials = 100
G.mat <- matrix(NA,  n_trials, n_gen)
for (i in 1:n_trials)
    {
        N = 1000
        A <- fwm(N, n_gen, mu = 4/N)
        G.vec <- apply(A, 1, get.G)
        G.mat[i,] <- G.vec
        }
matplot(t(G.mat), type = "l")

G.bar <- apply(G.mat, 2, mean)
lines(G.bar, lwd = 4)
## cool plot, why is it about .11?

1/(1 + 8)
## [1] 0.1111111

abline(h = 1/9, lty = 3, col = "yellow", lwd = 5)

## Gillespie tells us that Gbar is supposed to be 1 / (1 + 4*Ne*mu)
## How does 4*Ne*mu = 8?
## Well, we have K*mu = 4 and since K = 2*Ne, Ne = = K/2
## (maybe)

#######################################################
## FW babyname simulation of equilibrium frequencies ##
#######################################################

N = 1000 ## Not sure if this is (w)right :)
mu = 4/N ## [1] 0.004
theta = N*mu ## [1] 4 ## H&B's "best fit"


## Now we can simulate babyynmes
n_gen = 1001
N = 1000
## set.seed(1)
## A <- fwm2(N, n_gen, mu = 4/N)
A <- fwm2(N, n_gen, mu = 8/N)
## ok, lets do power law plot of this
x <- A[1001,]
tt <- table(x)
ptt <- prop.table(tt)
my.breaks <- c(0, 2^(0:11)/10000)
p <- ptt
## bar[, p := round(prop.table(N),5)]
##bar[, pcat := cut(p, breaks =  my.breaks, right = F, include.lowest = T)]
pcat = cut(p, breaks =  my.breaks, right = F, include.lowest = T)
out <- unclass(prop.table(table(bar$pcat)))
out.hat <- unclass(prop.table(table(pcat)))
my.x <- my.breaks[-length(my.breaks)] + diff(my.breaks)/2
plot(my.x, out, log = "xy")
lines(my.x, out.hat)




