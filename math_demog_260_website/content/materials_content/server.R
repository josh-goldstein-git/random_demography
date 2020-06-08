## server.R for exponential_growth_with_poisson_heterogeneity

library(shiny)

## functions used


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



exponential_growth_with_poisson_heterogeneity_fun <-
    function(r0 = .03,
             lambda= 2,
##             a_string = NULL, ## from radioButton
             a = .01,
             n_subpops = 11,
             t_years = 100,
             out_return = FALSE,
             as.pdf = FALSE)
{

##    if(!is.null(a_string))
##        a = as.numeric(a_string)


    ##  let r_i = r0 - s_i * a
    ## with s_i ~ poisson( mean = lambda)
    s.values = 0:(n_subpops - 1) ## e.g., 0 0:10
    f.of.s <- dpois(x = s.values, lambda)
    ## check for enough terms to cover mass of the Poisson
    my.sum <-  sum(f.of.s)
    if (my.sum < .99)
    {
        stop("Error: increase n_subpops for proper distribution summing close to 1.0")
    }

    r_i = r0 - s.values*a

    ## now simulate pop growth
    tt = 0:t_years
    Ki.mat <- matrix(NA, length(r_i), length(tt)) ## a matrix storing the size of each subpops
    for (i in 1:length(r_i))
    {
        Ki.mat[i,] <- f.of.s[i] * exp(r_i[i] * tt)
    }
    K.tot.vec <- apply(Ki.mat, 2, sum)
    K.tot.vec.norm <- K.tot.vec/K.tot.vec[1] ## normalize so that otal pop size starts at 1.0
    R.obs <- center.diff(log(K.tot.vec.norm))

    ## get variance of growth rates r_i
    Pi.mat <- prop.table(Ki.mat, 2) ## matrix of proportions in each subpop
    var.vec <- (rbind(r_i)^2 %*%  Pi.mat) - (rbind(r_i) %*%  Pi.mat)^2

    ## analytic solution
    K.tot.vec.hat <- exp(r0*tt) * exp(-lambda * (1 - exp(-a * tt)))
    R.hat <- r0 - a * lambda * exp(-a * tt)
    out <- list(Ki.mat = Ki.mat,
                K.tot.vec = K.tot.vec,
                K.tot.vec.norm = K.tot.vec.norm,
                K.tot.vec.hat = K.tot.vec.hat,
                R.hat = R.hat,
                var.vec = drop(var.vec))

    ## plot results
    if(as.pdf == TRUE)
    {
        pdf("kens_heterogeneity_approach_figure.pdf")
    }
    par(mfrow = c(2,2))
    plot(r_i, f.of.s, ## xlim = c(-.15, .15),
         ylab = "Share of world population",
         xlab = "Growth of sub-population, r_i",
         axes = F)
    title("Distribution of sub-population growth rates")
    axis(1, at = seq(-.15, .15, .01), cex.axis = .7, tick = TRUE)
##     axis(1, at = 0, labels = F, tick = TRUE)
    axis(2)
     abline(v = 0, col = "grey", lwd = 2)
    segments(x0 = r_i, x1 = r_i,
             y0 = 0, y1 = f.of.s,
             col = 1:6, lty = 1:6)
    text(x = r0, y = f.of.s[1], "r0", pos = 3)

    matplot(tt, t(Ki.mat), type = "l",
            ylab = "Size of sub-population")
    ## matplot(tt, t(Ki.mat), type = "l", log = "y")
    text(x = 80, y = Ki.mat[1,80], "r0", pos = 3)
    title("Sub-population sizes")

    plot(tt, K.tot.vec.norm, type = "l",
         ylab = "Total population size")
    title("Aggregate population size")
    lines(tt, K.tot.vec.norm, lty = 2, col = "red", lwd = 3)
    legend("topleft", legend = c("observed", "formula"),
           bty = "n",
           col = c("black", "red"), lty = c(1, 2), lwd = c(1,3))
### bingo. Everything works now!!
##     source("~/Documents/tempo_repos/optim_clean/utility_functions.R")


    cat(file = stderr(), "f.of.s", f.of.s)
    plot(tt, R.obs, type = "l",
         ylim = c(-.1, .1),
         ## min(r_i), max(r_i) * 2),
         ylab = "Growth rate")
    title("Aggregate growth rate")
##     R.hat <- r0 - a * lambda * exp(-a * tt)
    lines(tt, R.hat, lwd = 3, col = "red", lty =2)
    abline(h = r0)
    text(x = 0, y = r0, "r0", pos = 3)
##     legend("topleft", legend = c("observed", "formula"),
##           col = c("black", "red"), lty = c(1, 2), lwd = c(1,3))
    if(as.pdf == TRUE)
    {
        dev.off()
        system("open kens*figure.pdf") ## this may work only on apple mac (uses "open")
    }
    if(out_return == TRUE)
        return(out)
}


## use shorter name
poisson_het <-  exponential_growth_with_poisson_heterogeneity_fun

                                        # Define server logic required to draw a histogram
shinyServer(function(input, output) {

    ## Expression that generates our plot. The expression is
    ## wrapped in a call to renderPlot to indicate that:
    ##
    ##  1) It is "reactive" and therefore should re-execute automatically
    ##     when inputs change
    ##  2) Its output type is a plot
##     print(input$lambda)
    output$poissonPlot <- renderPlot({
         ## this is code for doing whole thing
        poisson_het(r0 = input$r0,
                    lambda = input$lambda,
##                    a_string = input$a_string,
                    n_subpops = input$n_subpops,
                    t_years = input$t_years)
    })
})
