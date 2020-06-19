## can use RCurl to download files, perhaps

get.year <- function(x) {as.numeric(names(x))}

logit <- function(p) { log (p  / (1-p) ) }
antilogit <- function(alpha) { 1 / (1 + exp(-alpha))}
last <- function(x){x[length(x)]}

my.cumsum <- function(x)
{
  ## drop leading NAs before cumsum
  ## cumsum will make everything NAs

  ## > my.cumsum(c(1,NA,1))
  ## [1]  1 NA NA
  ## > my.cumsum(c(NA,1,1))
  ## [1] NA  1  2
  ## > my.cumsum(c(NA,1,2, 3, NA, 1))
  ## [1] NA NA NA NA NA NA

    x.na.elements <- seq(x)[is.na(x)]
    x.nona.elements <- seq(x)[!is.na(x)]
    x.nona <- x[x.nona.elements]
    ## if leading NAs, then min(x.nona.elements)-1 = max(x.na.elements)
    if (min(x.nona.elements)-1 == max(x.na.elements))
    {
        out = x
        out[x.na.elements] <- NA
        out[x.nona.elements] <- cumsum(x.nona)
    }
    if (min(x.nona.elements)-1 != max(x.na.elements))
        out = cumsum(x)

    return(out)
 }


     ## check on location of NAs


get.moment <- function(fx, n.moment, x = NULL, n = NULL)
{
    ## compute moment from density fx
    ##    sum(((x+n/2)^n.moment)*fx)/sum(fx)

    ## usage:
    ## set.seed = 1
    ## y = rnorm(1000000, mean = 100, sd = 10)
    ## fy = table(floor(y))
    ## my.x2.bar = get.moment(fy, 2)


    ## assign age x if not given
    if (is.null(x))
    {
        ## if names, use
        if(!is.null(names(fx)))
            x = as.numeric(names(fx))
        ## if no names,  default 0, 1, ..
        if(is.null(names(fx)))
            x = seq(fx) - 1
    }
    ## assume open interval is same length as next to last
    if(is.null(n))
	n = c(diff(x), diff(x)[length(diff(x))])


    moment = sum(((x+n/2)^n.moment)*fx)/sum(fx)
    return(moment)
}

get.mean <- function(fx, x = NULL, n = NULL)
{
    get.moment(fx, n.moment = 1, x, n)
}
get.var <- function(fx, x = NULL, n = NULL)
{
    x.bar =  get.mean(fx, x, n)
    x2.bar = get.moment(fx, n.moment = 2, x, n)
    s2 = x2.bar - x.bar^2
    return(s2)
}
get.sd <- function(fx, x = NULL, n = NULL)
{
    s2 = get.var(fx, x = NULL, n = NULL)
    s = sqrt(s2)
    return(s)
}

## check to make sure these work

set.seed = 1
y = rnorm(1000000, mean = 100, sd = 10)
fy = table(floor(y))
get.mean(fy)
mean(y)
get.mean(fy)
var(y)
get.var(fy)
sd(y)
get.sd(fy)
## looks great


## try doing numerical derivatives using "filter"

if (0) {
x <- 1:10
x3 <- x^3
first.d <- filter(x3, -c(-1/2, 0, 1/2))
first.d.4 <- filter(x3,-c(1/12,-2/3,  0, 2/3, -1/12))
second.d <- filter(x3, c(1,-2,1))
second.d.4 <- filter(x3, c(-1/12, 4/3, -5/2, 4/3, -1/12))
}

second.diff <- function(x) {
    my.filter <- c(-1/12, 4/3, -5/2, 4/3, -1/12)
    second.d <- filter(x, my.filter)
    second.d
}
first.diff <- function(x) {
    my.filter <- -c(1/12,-2/3,  0, 2/3, -1/12)
    first.d <- filter(x, my.filter)
    first.d
}


center.diff <- function(x, end.fill = F)
{
    ## approximate derivatives with discrete data by taking central
    ## differences d(x) = ([f(x+1) - f(x)] + [f(x) - f(x-1)])/2
    ##                  = [f(x+1) - f(x-1)]/2
    ## if end.fill = T, then first and last differences are not centered.

    ## useful for Bongaarts-Feeney, Gompertz fertility, and other
    ## fitting of models that are developed in continuous time and
    ## involve derivatives

    forward.diff = c(diff(x), NA)
    backward.diff = c(NA, diff(x))
    out = (forward.diff + backward.diff)/2
    if(end.fill)
    {
        out[1] = diff(x)[1]
        out[length(out)] = diff(x)[length(diff(x))]
    }
    ## preserve names if exist
    if(!is.null(names(x)))
        names(out) = names(x)
    return(out)
}

remove.plusminus <- function(character.vector)
{
    ## replaces + and - with nothing
    ## e.g., "99+" becomes "99"
    out = gsub("[-+]", "", character.vector)
    return(out)
}



get.lexis.mat <- function(age.vec, time.vec, rate.vec)
{
    ## converts long to wide format, to create an "age" x "time"
    ## (cohort or period) matrix

    ## to make go faster we loop only over time, but we do check each
    ## time to make sure that all the ages are there and that they are
    ## in the right order

    age = names(table(age.vec))         # can leave as character
    time = names(table(time.vec))
    Fat = matrix(NA, length(age), length(time))

    dimnames(Fat) <- list(age, time)
    ## loop through each time period
    for (i in 1:length(time))
    {
        this.time = time[i]
        ## check to make sure all ages there and in right order
        tmp = time.vec == this.time
        all.there = identical(as.character(age.vec[tmp]), age)
        if (all.there)
        {
            Fat[,i] <- rate.vec[tmp]
        }
        if (!all.there)
            print(paste("Error for time", this.time))
    }
    return(Fat)
}

read.table.hfd <- function(filename, ...)
{
    df = read.table(file = filename, header = T,
    skip = 2, na.strings = ".", as.is = T, ...)
    return(df)
}



hfd.rate2array <- function(filename, var = "ASFR", code = NULL,
                           age.str = "Age")
{
    ## converts a "long" hfd rate file with info for multiple parities
    ## to an array using get.lexis.mat for each parity
    ## if not broken down by parity, still works, but returns matrix

    ## usage examples:
    ## BB = hfd.rate2array("../Fertility_Data/NLDasfrRRbo.txt")
    ## B = hfd.rate2array("../Fertility_Data/NLDasfrRR.txt")

    if (0) {
    code = "NLD"
    var = "ASFR"
    filename = "../hfd/asfrRRbo.txt"
}
    df = read.table.hfd(filename)
    if (!is.null(code))
    {
        ## then we have combined file with country codes,
        ## select only the one with the code
        df <- df[df$Code == code,]
    }


    df.names = names(df)

    ## figure out how many parities
    parity.str <- df.names[grep(var, df.names)]

    time.vec = df$Year
##    age.vec = df$Age
    age.vec = df[[age.str]]

    age = names(table(age.vec))
    time = names(table(time.vec))

    ## create array
    AAA <- array(NA, c(length(age), length(time), length(parity.str)))
    dimnames(AAA) = list(age, time, parity.str)

    ## loop through parities, filling array
    for (i in 1:length(parity.str))
    {
        this.parity.str = parity.str[i]
        this.rate.vec = df[[this.parity.str]]
        AAA[,,i] <- get.lexis.mat(age.vec, time.vec, this.rate.vec)
    }
    if (dim(AAA)[3] > 1)                # if broken down by birth order
        return(AAA)
    if (dim(AAA)[3] == 1)               # if no birth order return matrix
        return(AAA[,,1])
}




per2coh <- function(Aper)
{
  ## converts a period matrix into a cohort matrix
  ## Note: requires dimnames
  if(is.null(dimnames(Aper)))
    print("error in per2coh: dimnames of period matrix required")
  ##     Aper <- matrix(1:6, 2, 3)
  ##     dimnames(Aper) <- list(0:1, 0:2)
    per.names <- colnames(Aper)
    age.names <- rownames(Aper)
    periods <- as.numeric(per.names)
    ages <- as.numeric(age.names)

    cohorts <- (min(periods)-max(ages)):(max(periods)+max(ages)-1)

    Acoh <- matrix(NA, length(ages), length(cohorts))
    dimnames(Acoh) <- list(ages, cohorts)
    for (j in 1:length(periods))
        for (i in 1:length(ages))
        {
            this.per <- periods[j]
            this.age <- ages[i]
            this.coh = this.per - this.age
            this.coh.name <- as.character(this.coh)
            Acoh[i,this.coh.name] <- Aper[i,j]
        }
    return(Acoh)
}
all.na <- function(x)
{
    all(is.na(x))
}
remove.na.columns <- function(A)
{
    A[,!apply(A, 2, all.na)]
}
remove.na.rows <- function(A)
{
    A[!apply(A, 1, all.na),]
}

coh2per <- function(Acoh)
{
  ## converts a cohort matrix to a period matrix
  ## rm.na.cols is useful to perserve size of period matrix
  ## e.g.
  if(0) {
    Aper <- matrix(1:6, 2, 3)
    dimnames(Aper) <- list(0:1, 0:2)
    Acoh <- per2coh(Aper)
    Aper.recreate <- remove.na.columns(coh2per(Acoh))
  }

  ## Note: requires dimnames

  coh.names <- colnames(Acoh)
  age.names <- rownames(Acoh)
  cohorts <- as.numeric(coh.names)
  ages <- as.numeric(age.names)

  if(is.null(dimnames(Acoh)))
    print("error in coh2per: dimnames of cohort matrix required")

  ##     periods <- (-length(ages) + 1):(max(cohorts) + length(ages)-1)
  ##    periods <- min(cohorts):(max(cohorts) + max(ages) - 1)
  periods <- min(cohorts):(max(cohorts) + max(ages))
  ## cohorts <- (min(periods)-max(ages)):(max(periods)+max(ages)-1)

  Aper <- matrix(NA, length(ages), length(periods))
  dimnames(Aper) <- list(ages, periods)
  ##  print(dimnames(Aper))
  for (j in 1:length(cohorts))
    for (i in 1:length(ages))
      {

        this.coh <- cohorts[j]
        this.age <- ages[i]
        this.per = this.coh + this.age
        this.per.name <- as.character(this.per)
        Aper[i,this.per.name] <- Acoh[i,j]
      }
  return(Aper)
}


convert.parity.array <- function(AAA)
{
    ## add together the parities bigger than 3 to create a 3+ category
    parities <- dimnames(AAA)[[3]]
    parities.3plus <- parities %in% c("ASFR3", "ASFR4", "ASFR5p")
    A.3plus <- apply(AAA[,,parities.3plus], 1:2, sum)
    new.parities <- c("ASFR", "ASFR1", "ASFR2", "ASFR3p")
    AAA.out <- array(NA, c(nrow(AAA), ncol(AAA),
                           length(new.parities)))
    dimnames(AAA.out) <- list(rownames(AAA), colnames(AAA),
                              new.parities)
    AAA.out[,,1:3] <- AAA[,,1:3]
    AAA.out[,,"ASFR3p"] <- A.3plus
    if (sum(AAA.out) != sum(AAA))
        print("error in convert.parity.array: \n sum(AAA.out) != sum(AAA))")
    return(AAA.out)
}

