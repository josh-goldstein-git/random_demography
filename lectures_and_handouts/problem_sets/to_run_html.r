# To compile into html the relevant files are:
# - _bookdown.yml: change name of book_filename
# - _output.yml: change the settings or add them to the YAML in the markdown.

# Code to render as html
# gitbook(fig_caption = TRUE, number_sections = TRUE,
#         self_contained = FALSE, lib_dir = "libs",
#         pandoc_args = NULL,  template = "default",
#         split_by =  "none",
#         split_bib = TRUE, config = list(), table_css = TRUE)
#
# render_book("index.Rmd", "bookdown::gitbook")

render_book("index.Rmd","bookdown::gitbook",
            clean = TRUE, envir = parent.frame(),
            clean_envir = !interactive(), output_dir = NULL,
            new_session = NA, preview = FALSE,
            config_file = "_bookdown.yml")

# html_chapters("index.Rmd",toc = TRUE, number_sections = TRUE,
#               fig_caption = TRUE, lib_dir = "libs",
#               template = "/hdir/0/andreamg/R/x86_64-pc-linux-gnu-library/3.6/bookdown/templates/default.html",
#               pandoc_args = NULL,
#               base_format = rmarkdown::html_document,
#               split_bib = TRUE, page_builder = build_chapter,
#               split_by =  "none")
