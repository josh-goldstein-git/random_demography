hgomp <- function(x, b, a)
{
  a * exp(b * x)
}
pgomp <- function(x,b, a)
{
  Hx <- (a/b) * (exp(b*x) - 1)
  lx <- exp(-Hx)
  1-lx
}
dgomp <- function(x, b, a)
{
  hx <- hgomp(x,b,a)
  lx <- 1-pgomp(x,b,a)
  dx <- hx*lx
  dx
}
idfgomp <- function(u, b, a)
{
  x = (1/b) * log(-log(1-u)*b/a + 1)
  x
}
## idfgomp(u = .5, b = .1, a = .001)
rgomp <- function(N, b, a)
{
  u <- runif(N)
  x <- idfgomp(u, b, a)
  x
}

get.Dx <- function(y)
{
  ## counts number of x in single year age groups
  ## including zeros when there's no one
  ## (note: built-in "table()" won't do this :(
  x <- 0:max(floor(y))
  y.fac <- factor(floor(y), levels = x)
  Dx <- tabulate(y.fac)
  names(Dx) <- x
  return(Dx)
}
