#Frailty simulation function using gompertz mortality
frailty_sim <- function( N, z, base.a, base.b ){
  ## (1) simulate ages at death from h0*z  using the gompertz as our baseline##
  ## note: we call the continuous ages of death "y"
  ## but we'll make a table of deaths at age "x" and
  ## using life table notation call the count "Dx"
  y <- rgomp(N,
             b = base.b, ## doesn't vary
             a = base.a * z) ## multiplicative fixed frailty

  ## (2) Lifetables: first define age at death as floor(y) and then
  ## make a table of deaths at each age ("Dx")

  Dx <- get.Dx(y)
  x <- as.numeric(names(Dx))
  lx <- rev(cumsum(rev(Dx))) ## lx by reverse-survival
  lxpn <- c(lx[-1], 0) ## Person-years as average of adjacent lx
  Lx <- (lx + lxpn)/2
  mx <- Dx/Lx ## Hazards
  Tx <- rev(cumsum(rev(Lx))) ## Remaining person-years
  ex <- Tx/lx ## Life expectancy at age x

  ## Baseline lifetable
  lx.base <- N * (1- pgomp(x, b = base.b, a = base.a))
  Dx.base <- round(-diff(c(lx.base,0)))
  mx.base <- base.a * exp(base.b * (x + .5)) ## x + .5
  lxpn.base <- c(lx.base[-1], 0)
  Lx.base <- (lx.base + lxpn.base)/2
  Tx.base <- rev(cumsum(rev(Lx.base)))
  ex.base <- Tx.base/lx.base


  # exported tables
  lifetables <- list()
  lifetables$sim <- y
  lifetables$z <- z
  lifetables$baseline <- tibble(Dx.base, lx.base,lxpn.base, Lx.base, mx.base, Tx.base, ex.base)
  lifetables$frailty <- tibble(x,Dx, lx,lxpn, Lx, mx, Tx, ex)
  return(lifetables)

}
