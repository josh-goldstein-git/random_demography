## Tempo adjustment of US fertility using HMD data and
## Bongaarts-Feeney formula

library(data.table)
source("~/Documents/tempo_repos/optim_clean/tempo_functions.R")
source("~/Documents/tempo_repos/optim_clean/utility_functions.R")


## (1) read in data and format into an array
## (2) fit bongaarts feeney without birth order
## (3) fit bongaarts feeney with birth order

## (1) read in data and format into an array

## age specific fertility rates by birth order for all countries and times
## RR means "rectangles" on Lexis surface
dt <- fread("~/Documents/hfd/Files/zip_w/asfrRRbo.txt")
dt <- dt[Code == "USA"] ## keep only US
## keep only ages 15 to 49
dt <- dt[Age %in% 15:49]
print(dt)

## put all order fertiility into a matrix

fat <- dt[, xtabs(ASFR ~ Age + Year)]
fat <- as.matrix(unclass(fat)) ## may cause problems, we'll see

fat1 <- dt[, xtabs(ASFR1 ~ Age + Year)]
fat2 <- dt[, xtabs(ASFR2 ~ Age + Year)]
fat3 <- dt[, xtabs(ASFR3 ~ Age + Year)]
fat4 <- dt[, xtabs(ASFR4 ~ Age + Year)]
fat5p <- dt[, xtabs(ASFR5p ~ Age + Year)]

year.vec <- colnames(fat)
age.vec <- rownames(fat)
parity.vec <- c("all", 1:5)
fat.array <- array(NA, dim = c(nrow(fat), ncol(fat), length(parity.vec)))
dimnames(fat.array) <- list(age.vec, year.vec, parity.vec)
fat.array[,,"all"] <- fat
fat.array[,,"1"] <- fat1
fat.array[,,"2"] <- fat2
fat.array[,,"3"] <- fat3
fat.array[,,"4"] <- fat4
fat.array[,,"5"] <- fat5p

## (2) fit bongaarts feeney

tfr.vec <- colSums(fat)

## (2a) by hand
mu.vec <- apply(fat, 2, get.mean)
rt.vec <- center.diff(mu.vec)
adj.tfr.vec <- tfr.vec / (1 - rt.vec)

par(mfrow = c(3,1))
plot(year, mu.vec)
plot(year, rt.vec)
abline(h =0)
plot(year.vec, tfr.vec, type = "l")
lines(year.vec, adj.tfr.vec, lty = 2)
abline(v = c(1945, 2008))
## we see fertility since 1980 has been depressed by postponment
## we see weird dynamics around end of WW2 and great recession.
## what's going on?

par(mfrow = c(1,1))
plot(year.vec, tfr.vec, type = "l")
lines(year.vec, adj.tfr.vec, lty = 2)
abline(v = c(1945, 2008))
abline(h = seq(1.5, 2.5, .1), col = "grey", lty = 3)

## can also use function to fit

adj.tfr.vec.from.fun <- bf.fit.simple(fat)$tfr.star
lines(year.vec, adj.tfr.vec.from.fun, col = "red")

## Now let's look at turbulence around WW2
par(mfrow = c(2,2))
plot(age.vec, fat[,"1944"], type = "l", ylim = c(0, .23),
     ylab = "f(a)",
     xlab = "age a"
     )
lines(age.vec, fat[,"1945"], type = "l", col = "red")
lines(age.vec, fat[,"1946"], type = "l", col = 'orange')
lines(age.vec, fat[,"1947"], type = "l", col = "blue")
legend("topright",
       legend = 1944:1947,
       col = c("black", "red", "orange", "blue"),
       lty = 1)
title("Age specific fertility")
##
plot(1943:1947, mu.vec[paste(1943:1947)],
     ylab = "mu(t)",
     xlab = "t",
     col = c("black", "black", "red", "orange", "blue"),
     pch = 19)
title("Mean ages")
##
plot(1943:1947, rt.vec[paste(1943:1947)],
     ylab = "r(t)",
     xlab = "t",
     col = c("black", "black", "red", "orange", "blue"),
     pch = 19)
title("Changes in mean, centered")
##
plot(1943:1947, tfr.vec[paste(1943:1947)],
     ylab = "tfr",
     xlab = "t",
     ylim = c(1, 4),
     type = "l")
lines(1943:1947, adj.tfr.vec[paste(1943:1947)],
      lty = 2)
title("TFR and adjTFR")

## From 1945 to 1946, fertility goes up a lot, but more at younger
## ages. So mean goes down. BF adjustment over-compensates, and has
## quantum declining.

## What's happening from 1944-45?


## (3) fit bongaarts feeney with birth order

out <- bf.fit(fat.array)

adj.tfr.bo.vec <-  out$tfr.star.out[, "bf.tfr.star"]

## pdf("usa_tempo_fig.pdf")
par(mfrow = c(1,1))
plot(year.vec, tfr.vec, type = "l", lwd = 2)
lines(year.vec, adj.tfr.vec, lty = 2)
lines(year.vec, adj.tfr.bo.vec, lty = 1, lwd = 2, col = "red")
## let's check against hfd
hfd.adj.dt <- fread("~/Documents/hfd/Files/zip_w/adjtfrRR.txt",
                skip = 2)
hfd.adj.dt <- hfd.adj.dt[Code == "USA"]
hfd.adj.dt[, points(Year, adjTFR, col = "red", cex = .8)]
legend("topright",
       c("tfr", "tfr* (all parities)", "tfr* (by parity)"),
       col = c("black", "black", "red"),
       lty = c(1, 2, 1),
       lwd = c(2,1,2))
hfd.adj.dt[, lines(Year, filter(adjTFR, rep(1/7, 7)), col = "red", lwd = 4)]

## dev.off()

## Q. Does taking birth order into account smooth WW2 turbulence?
## Q. Does taking birth order into account smooth Recession turbulence?
## Q. Does taking birth order into account, retell the baby boom?
## Q. Does taking birth order into account, retell the baby bust?

