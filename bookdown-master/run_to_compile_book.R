library(bookdown)
<<<<<<< HEAD
library(data.table)
setwd("/hdir/0/fmenares/Book/bookdown-master/")
=======
# setwd("/hdir/0/fmenares/Book/bookdown-master/")
setwd("/hdir/0/andreamg/Year2_2019_2020/Random_demography/random_demography/bookdown-master/")

>>>>>>> f0614114de38ec3c060c45ecec38f3138c202de1
render_book("index.Rmd","bookdown::gitbook",
            clean = TRUE, envir = parent.frame(),
            clean_envir = !interactive(), output_dir = NULL,
            new_session = NA, preview = FALSE,
            config_file = "_bookdown.yml")



