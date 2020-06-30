render_book("index.Rmd","bookdown::gitbook",
            clean = TRUE, envir = parent.frame(),
            clean_envir = !interactive(), output_dir = NULL,
            new_session = NA, preview = FALSE,
            config_file = "_bookdown.yml")