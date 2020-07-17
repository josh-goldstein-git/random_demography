library(bookdown)
# setwd("/hdir/0/fmenares/Book/bookdown-master/")
setwd("/hdir/0/andreamg/Year2_2019_2020/Random_demography/random_demography/bookdown-master/")

render_book("index.Rmd","bookdown::gitbook",
            clean = TRUE, envir = parent.frame(),
            clean_envir = !interactive(), output_dir = NULL,
            new_session = NA, preview = FALSE,
            config_file = "_bookdown.yml")


