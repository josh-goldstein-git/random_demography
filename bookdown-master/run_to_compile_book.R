library(bookdown)
#setwd("/hdir/0/fmenares/Book/bookdown-master/")
setwd("/hdir/0/andreamg/Year2_2019_2020/Random_demography/random_demography/bookdown-master/")
render_book("index.Rmd","bookdown::gitbook",
            clean = TRUE, envir = parent.frame(),
            clean_envir = !interactive(), output_dir = NULL,
            new_session = NA, preview = FALSE,
            config_file = "_bookdown.yml")
#In order to run this file you have to edit the directory (setwd) and the output dir
#in the _bookdown.yml file

#In order to modify the chapters included, you must edit _bookdown.yml list. 

# To preview a chapter just run:
preview_chapter('09-coalescent.Rmd')
