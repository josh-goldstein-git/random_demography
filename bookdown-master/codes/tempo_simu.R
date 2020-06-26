## Shift simulator

## We simulate births according to their original timing and then
## shift them by age and time according to a continuous shift function.

## We then show that changes in the mean age of the shifted births can
## be used to recover the original birth counts.

## Finally, we consider the case when the population consists of two
## ("hidden") sub-populations. The first population is identical to
## the population above. The second population are "early"
## births. Here we assume there is no change in timing, but there is a
## decline in the number of early births.

## We see what happens to the mean ages in this situation and whether
## we can recover the orignal births.


library(data.table)
source("/hdir/0/fmenares/Book/bookdown-master/codes/utility_functions.R")
source("/hdir/0/fmenares/Book/bookdown-master/codes//tempo_functions.R")


million = 10^6
thousand = 10^3
N <-  1 * million
year.vec <- 1990:2020

#######################################
## simulate originally timed births  ##
#######################################

## we'll assume normal with age and uniform with time
x <- rnorm(N, mean = 30, sd = 4)
t <- runif(N, min = min(year.vec), max(year.vec)+1)
dt <- data.table(x, t)

par(mfrow = c(1,1))
if (N == 1000)
{
    dt[, plot(t, x, cex = .6)]
    }


#########################
## shifting the births ##
#########################

## we'll assume a continuous S-shaped cumulative shift
## shift.t <- 2*plogis(seq(-5, 5, length.out = length(year.vec)))
 shift.t <- sin((year.vec - 1990)/2.5)


plot(year.vec, shift.t, type = "o",
     main = "Cumulative shifts")



## include shifts as a column
Rt <- approx(x = year.vec, y = shift.t, xout = t)$y
dt[, Rt := approx(x = year.vec, y = shift.t, xout = t)$y]
## calculate shifted times and ages of births
dt[, t.obs := t + Rt]
dt[, x.obs := x + Rt]
## retain only the original time window (for convenience)
dt <- dt[floor(t.obs) %in% year.vec]

if (N == 1000)
{
    dt[, plot(t, x, cex = .6, col = "grey")]
    dt[, points(t.obs, x.obs, cex = .6, col = "orange", pch = 19)]
    }

##########################################
## observed births counts and mean ages ##
##########################################

out <- dt[, .(Bt = .N,                   ## count the events
              mut = mean(x.obs)),      ## get mean age
          by = .(year = floor(t.obs))] ## by whole years
out <- out[order(year)]

############################################
## change in mean age and adjusted counts ##
############################################

out[, rt.hat := center.diff(mut, end.fill = T)]
out[, Rt.hat := cumsum(rt.hat)]
out[, Bt.adj := Bt / (1 - rt.hat)]


######################
## plot the results ##
######################

par(mfrow = c(2,2))
out[, plot(year, Bt, ylim = c(.8, 1.2) * range(Bt),
           main = "Observed Births")]
out[, plot(year, mut,
           main = "Mean age of birth")]
out[, plot(year, center.diff(mut),
           main = "Change in mean age of birth")]
abline(h = 0)
## observed, adjusted, and original births
Bt.orig.vec <- dt[, table(floor(t))]
out[, plot(year, Bt, ylim = c(.8, 1.5) * range(Bt),
           main = "Observed and Adjusted Births")]
out[, lines(year, Bt.adj, col = "red")]
points(names(Bt.orig.vec), Bt.orig.vec, col = "grey")
legend("top", c("observed", "adjusted", "original"),
       pch = c(1,-1,1), lty = c(-1, 1,-1),
       col = c("black", "red", "grey"),
       bty = "n")

## function version of tabulating and plotting

tempo_simu_plot_fun <- function(dt)
{

    ## requires x.obs and t.obs and
    ## (optionally) t, the original unshifted birth times

##########################################
## observed births counts and mean ages ##
##########################################

out <- dt[, .(Bt = .N,                   ## count the events
              mut = mean(x.obs)),      ## get mean age
          by = .(year = floor(t.obs))] ## by whole years
out <- out[order(year)]

############################################
## change in mean age and adjusted counts ##
############################################

out[, rt.hat := center.diff(mut, end.fill = T)]
out[, Rt.hat := cumsum(rt.hat)]
out[, Bt.adj := Bt / (1 - rt.hat)]


######################
## plot the results ##
######################

par(mfrow = c(2,2))
out[, plot(year, Bt, ylim = c(.8, 1.2) * range(Bt),
           main = "Observed Births")]
out[, plot(year, mut,
           main = "Mean age of birth")]
out[, plot(year, center.diff(mut),
           main = "Change in mean age of birth")]
## observed, adjusted, and original births
Bt.orig.vec <- dt[, table(floor(t))]
out[, plot(year, Bt, ylim = c(.8, 1.5) * range(Bt),
           main = "Observed and Adjusted Births")]
out[, lines(year, Bt.adj, col = "red")]
points(names(Bt.orig.vec), Bt.orig.vec, col = "grey")
legend("top", c("observed", "adjusted", "original"),
       pch = c(1,-1,1), lty = c(-1, 1,-1),
       col = c("black", "red", "grey"),
       bty = "n")
}

tempo_simu_plot_fun(dt)

## In-class exercises:

## (1) try with N of 4 million -- does it still work? What happens?

## (2) try with a shift function that goes up and down

 shift.t <- sin((year.vec - 1990)/2.5)
 plot(year.vec, shift.t, type = "o")

## Are the adjusted counts ever LESS than the observed counts? If so, when?

## (3) If the cumulative shift was Rt = a + 0.1*t, what would be a
## formula for tempo-adjusted counts of births? Sketch the 4 panels
## without the computer and then check to see if you're right.


