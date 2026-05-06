library(shiny)
library(shinyFiles)
library(bslib)
library(DT)

app_dir <- normalizePath(getwd())

source(file.path(app_dir, "R_util", "conversion_funcs.R"))
source(file.path(app_dir, "R_util", "util_funcs.R"))
source(file.path(app_dir, "R_util", "genesets.R"))
source(file.path(app_dir, "R_util", "reading_funcs.R"))
source(file.path(app_dir, "R_analysis", "DEGSEA.R"))
source(file.path(app_dir, "R_analysis", "signature_projection.R"))
source(file.path(app_dir, "R_analysis", "anchored_GSEA.R"))
source(file.path(app_dir, "R_analysis", "OUTRIDER.R"))
source(file.path(app_dir, "R_analysis", "univariate_Cox.R"))
source(file.path(app_dir, "shiny_modules", "tutorials", "DEGSEA_tutorial.R"))
source(file.path(app_dir, "shiny_modules", "tutorials", "OUTRIDER_tutorial.R"))
source(file.path(app_dir, "shiny_modules", "mod_DEGSEA.R"))
source(file.path(app_dir, "shiny_modules", "mod_signature_projection.R"))
source(file.path(app_dir, "shiny_modules", "mod_anchored_GSEA.R"))
source(file.path(app_dir, "shiny_modules", "mod_OUTRIDER.R"))
source(file.path(app_dir, "shiny_modules", "mod_univariate_Cox.R"))

options(shiny.maxRequestSize = 4096 * 1024^2) # Max supported input file size

ui <- page_navbar(
  nav_panel("Welcome "),
  nav_panel("DEGSEA", mod_degsea_ui("degsea")),
  nav_panel("Signature Projection", mod_signature_proj_ui("sign_proj")),
  nav_panel("Anchored GSEA", mod_anchored_gsea_ui("anchored_gsea")),
  nav_panel("OUTRIDER", mod_outrider_ui("outrider")),
  nav_panel("Univariate Cox", mod_univariate_cox_ui("univariate_cox"))

#   nav_panel("GSVA", mod_gsva_ui("gsva"))  # later
)

server <- function(input, output, session) {
  # roots <- c(home = "~")
  roots <- c(home = Sys.getenv("SHINY_ROOT_PATH", "/browse"))
  shinyDirChoose(input, "dir", roots = roots)

  mod_degsea_server("degsea", roots = roots)
  mod_signature_proj_server("sign_proj", roots = roots)
  mod_anchored_gsea_server("anchored_gsea", roots = roots)
  mod_outrider_server("outrider", roots = roots)
  mod_univariate_cox_server("univariate_cox", roots = roots)

  # mod_gsva_server("gsva", roots = roots)
}

shinyApp(ui, server)