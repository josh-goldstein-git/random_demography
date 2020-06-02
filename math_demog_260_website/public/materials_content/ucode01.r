#    ucode01.R from  kcode01.R  R-Code for Hazard Plateaus Session of Demography 250 
#             20 February 2020;  Kenneth W. Wachter
#             Based on adir.barbi/gamgom.r 

####################################################################
######      Calculate the aggregate population hazard function predicted from 
######          Gompertz hazards for individuals combined with
######          a Gamma Distribution for fixed proportional frailty
######          with parameters matched to the Italian cohort data of
######          Barbi et al. (2018).
####################################################################

#      Parameters for Gompertz and Gamma 
    beta     <-  0.088            #  Gompertz slope parameter,  Barbi et al.
    alpha60  <-  0.01340671       #  Gompertz intercept parameter at age 60 
    asymp    <-  0.620            #  Observed asymptote, for comparisons in plots
    khat     <-  7.045455         #  Gamma shape parameter, to fit plateau 
    theta    <-  1/khat           #  Gamma rate parameter, for unit mean at 60

#      Set age range, starting to observe cohorts at age 60, reaching up to 130 
    xxx      <-  c(0:70)          #  set of values of x, years past 60
    age      <-  60  +  xxx       #  age  


###>>>     Insert the Gompertz formula to calculate individual hazards hx at x:

###  hx       <-  ???   

###>>>     Insert the formula to calculate individual cumulative hazard  Hx  

##    Hx       <-  ???


###>>>     Insert the formula for the aggregate population hazard mu_x with Gamma frailty:

##    mux       <-  ???

####>>>     Plot  the logarithm of the individual hazard hx versus age.
####>>>         Put in line for logarithm of the predicted population hazard  mux
####>>>         Put in line segment between x = 105 and x = 117 for
####>>>               the logarithm of the estimated plateau at  0.620  from Barbi et al.

####>>>     Does   mux approach its own plateau?
####>>>     Does the level of the plateau correctly match the observed level with these parameters?









#########################################################################
######     How does the predicted plateau from a Gamma Gompertz depend
######          on the shape parameter  k   and the Gompertz slope beta?
########################################################################


###>>>     Halve the value of k  and recompute  mux.
###>>>            What is the new plateau?

###>>>     Keeping this half-value of k,  halve the value of beta and recompute mux
###>>>            What is the new plateau?  

###>>>     Based on these or other tries,  propose a formula for the level of the
###>>>            plateau as a function of  k  and beta in the Gamma Gompertz model.









############################################################################
###   Plot the mean frailty in the original Gamma-Gompertz fit as it changes over age 
############################################################################

####>>>>    The mean frailty among survivors is the ratio of  mux to hx.    Why? 

####>>>>    Plot  mean frailty among survivors  as a function of age.

       
##      pdf("fig_meanfrail")
    plottitle <- "Mean Frailty from Gamma-Gompertz Model" 
    header <- paste("KWW " , date(), sep="   ")

    plot( age, mux/hx, xlim = c(60, 130), type = "l"  , xlab = "", ylab = "Mean Frailty" )
    abline( v = 105,  lty = 3     )                     # vertical age 105   

##       dev.off()

####>    What is the approximate mean frailty of survivors to 120?  

####>    Consider a woman whose individual frailty is equal to this value.
####>        When she was as young as 60,  what was her hazard?

     round( hx[1:10]/5, 3 )
 #            0.003 0.003 0.003 0.003 0.004 0.004 0.005 0.005 0.005 0.006


    plot( age[1:20], hx[1:20], type = "l" , ylim = c(0,0.10) )
    lines( age[1:20], (1/5)*hx[1:20]  )


############################################################################
###        Explore an example of a genetic process which can generate a 
###            plateau in individual hazard functions.  
###            It is based on the genetic evolutionary process of mutation 
###            accumulation.  
###        Begin with a case in which a deleterious mutant allele raises
###            the hazard function of those who carry it with onset
###            only in a single age interval.
###############################################################################

#                   Set parameter values 
    lambda   <-   0.080           #  background extrinsic hazard 
    delta    <-   0.002           #  increment to hazard from mutant allele
    epsilon  <-   0.0005          #  fixed cost 
    nu       <-   0.200/50        #  mutation rate per site per generation

#                   Define the range of ages  
    xxx      <-  c(0:50)          #  set of values of x, years past 20, up to 70
    age      <-  20  +  xxx       #  age  

#                   Find the fertility level that makes baseline NRR equal to 1.
    lx       <-  exp(-lambda*xxx) #  survivorship under background hazard 
    fert     <-  1/sum(lx)

#                   Set up a loop over choices of an age interval a at which
#                        the mutant allele affects the hazard.

    bigrho   <-  NULL             #  vector for holding values of rho  

    for  ( a in seq(xxx))      {

        hxa    <- lambda + 0*xxx + 0*epsilon        #  baseline hazard (no fixed cost)  
        hxa[a] <- hxa[a]  + delta                   #  add hazard increment in age interval a 

        Hxa   <- c(0, cumsum(hxa))[seq(hxa)]
        lxa   <- exp(-Hxa) 

####>>>>     Insert a formula for selective cost  Sa
####>>>>     Sa    <-  ???

####>>>>     Insert a formula for rho(a)  in terms of cost Sa and mutation rate nu
####>>>>     rho   <-  ????

        Sa    <- 1 - sum( fert*lxa)   
        rho   <- nu/Sa

        bigrho <- c(bigrho, rho)          #  keep the rho values in the vector bigrho 
                                   }      #  end of loop of choices of onset age a. 


####>>>>    Plot  bigrho as function of age 
####>>>>    Plot  log(bigrho)  as function of age
####>>>     Is there a stretch of ages that resemble an Gompertz exponential curve?

        plot(  age,  bigrho  )
        plot(  age,  log(bigrho)  )

#              Suppose an individual carries a sample of mutations of all types a,
#                   with mean proportional to the count rho(a) in the population.  
#                   The hazard on average would look like  hxx:

        hxx  <-  lambda +  bigrho*delta  

        plot ( age, hxx,  ylim = c(0,0.500) )     

####>>>>         Does such a hazard look like a Gompertz hazard?   A  Makeham hazard?


############################################################################
###        Proceed to a case in which a deleterious mutant allele imposes 
###            a fixed cost of size epsilon as well as adding an increment
###            of size delta to the hazard in age interval  a.
###############################################################################


#                 Repeat the loop over choices of an age interval a at which
#                        the mutant allele affects the hazard, this time with fixed cost..

    bigrho   <-  NULL             #  vector for holding values of rho  
    for  ( a in seq(xxx))      {

        hxa    <- lambda + 0*xxx + 0*epsilon        #  baseline hazard plus fixed cost epsilon
        hxa[a] <- hxa[a]  + delta                   #  add hazard increment in age interval a 

        Hxa   <- c(0, cumsum(hxa))[seq(hxa)]
        lxa   <- exp(-Hxa) 
        
        Sa    <- 1 - sum( fert*lxa)       #   Selective cost including allele at age a 
        rho   <- nu/Sa                    #   Equilibrium count of mutation a 

        bigrho <- c(bigrho, rho)          #  keep the rho values in the vector bigrho 
                                   }      #  end of loop of choices of onset age a. 


####>>>>    Now plot bigrho as function of age 
####>>>>       and plot log(bigrho)  as function of age
####>>>     Is there the appearance of a plateau?  

        plot(  age,  bigrho  )
        plot(  age,  log(bigrho)  )


#################################################################################
##################################################################################

  




 




