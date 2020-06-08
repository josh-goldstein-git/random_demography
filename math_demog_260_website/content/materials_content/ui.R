## ui.R for exponential_growth_with_poisson_heterogeneity

library(shiny)

shinyUI(fluidPage(
  # Application title
  titlePanel("Exponential Growth with Poisson heterogeneity"),
  h4("Choose levels and slide time"),
  p("The simulation starts in equilibrium.
To see the effect of exogenous changes in the starting population N, the wage function w(N), the birth function b(w), or the death function d(w), move the slider and hit the 'play' button the time slider."),
  p("For example, move the N slider to '8', hit 'play',  and watch the economy and population return to the steady state wage level with population '5'"),
  p("To replay more slowly, click on circular slider itself and use right arrow to advance (or left arrow to go back in time)."),
  # Sidebar with a slider input for the number of bins
  sidebarLayout(
    sidebarPanel(
      sliderInput("r0",
                  "Largest sub-population growth rate",
                  min = -.02,
                  max = .1,
                  value = .03),
      sliderInput("lambda",
                  "location parameter: r0 - lambda * a is mean growth rate",
                  min = 1,
                  max = 5,
                  value = 2),
      ## radioButtons("a_string",
      ##              "step size between sub-groups",
      ##              inline = T,
      ##              selected = "default",
      ##              c("0.001" = "tiny",
      ##                "0.01" = "default",
      ##                "0.1" = "huge")),
      sliderInput("n_subpops",
                  "Number of sub-groups",
                  min = 1,
                  max = 50,
                  value = 11),
      sliderInput("t_years",
                  "Number of years",
                  min = 10,
                  max = 500,
                  value = 100)
     ),

    ## Show a plot of the generated distribution
    mainPanel(
      plotOutput("poissonPlot"),
      h6("Note: new  conditions start at time 0")
      #,textOutput("Eq.b")
    )
  ) ## sidebarlayout(
) ## shinyUI(fluidpage(
) ## shinyUI(
