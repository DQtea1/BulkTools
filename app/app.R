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
source(file.path(app_dir, "R_analysis", "tximport_merge.R"))
source(file.path(app_dir, "R_analysis", "signatures_comparison.R"))
source(file.path(app_dir, "R_analysis", "ICA.R"))
source(file.path(app_dir, "shiny_modules", "tutorials", "DEGSEA_tutorial.R"))
source(file.path(app_dir, "shiny_modules", "tutorials", "OUTRIDER_tutorial.R"))
source(file.path(app_dir, "shiny_modules", "tutorials", "Signatures_comparison_tutorial.R"))
source(file.path(app_dir, "shiny_modules", "tutorials", "signature_projection_tutorial.R"))
source(file.path(app_dir, "shiny_modules", "tutorials", "ICA_tutorial.R"))
source(file.path(app_dir, "shiny_modules", "mod_Tximport.R"))
source(file.path(app_dir, "shiny_modules", "mod_DEGSEA.R"))
source(file.path(app_dir, "shiny_modules", "mod_signature_projection.R"))
source(file.path(app_dir, "shiny_modules", "mod_Signatures_comparison.R"))
source(file.path(app_dir, "shiny_modules", "mod_anchored_GSEA.R"))
source(file.path(app_dir, "shiny_modules", "mod_OUTRIDER.R"))
source(file.path(app_dir, "shiny_modules", "mod_univariate_Cox.R"))
source(file.path(app_dir, "shiny_modules", "mod_ICA.R"))

options(shiny.maxRequestSize = 4096 * 1024^2) # Max supported input file size

ui <- tagList(
  tags$head(tags$style(HTML("
    .responsive-plot-frame {
      height: clamp(320px, 72vh, 900px);
      width: 100%;
    }
    .responsive-plot-frame .shiny-plot-output,
    .responsive-plot-frame .shiny-image-output {
      width: 100% !important;
      height: 100% !important;
    }
    .responsive-plot-frame .shiny-image-output {
      display: flex;
      align-items: center;
      justify-content: center;
    }
    .responsive-plot-frame .shiny-image-output img {
      width: 100%;
      height: 100%;
      object-fit: contain;
    }
  "))),
  page_navbar(
    nav_panel("Welcome "),
    nav_panel("DEGSEA", mod_degsea_ui("degsea")),
    nav_panel("Signature Projection", mod_signature_proj_ui("sign_proj")),
    nav_panel("Anchored GSEA", mod_anchored_gsea_ui("anchored_gsea")),
    nav_panel("OUTRIDER", mod_outrider_ui("outrider")),
    nav_panel("Univariate Cox", mod_univariate_cox_ui("univariate_cox")),
    nav_panel("Tximport merge", mod_tximport_ui("tximport")),
    nav_panel("Signatures comparison", mod_signatures_comparison_ui("sig_comparison")),
    nav_panel("ICA", mod_ica_ui("ica"))

  #   nav_panel("GSVA", mod_gsva_ui("gsva"))  # later
  )
)

# Build the shinyFiles roots (the drives/folders the file & directory pickers
# can browse). The app runs in a container, so each root must be a path that is
# bind-mounted into the container (see the launchers). Roots are, in order:
#   - the primary browse root (host home by default)       -> SHINY_ROOT_PATH
#   - the output mount /out (if present)
#   - every directory mounted under /mounts/<label>        (external drives)
#   - explicit extras via SHINY_EXTRA_ROOTS="name=/path,name2=/path2"
# Falls back to "/" so nothing mounted is ever unreachable.
build_shiny_roots <- function() {
  roots <- character(0)

  add_root <- function(name, path) {
    if (is.null(path) || !nzchar(path) || !dir.exists(path)) return(invisible())
    name <- gsub("[^A-Za-z0-9_]+", "_", name)          # safe key for shinyFiles
    if (!nzchar(name)) name <- "root"
    base_name <- name; i <- 1L
    while (name %in% names(roots)) { i <- i + 1L; name <- paste0(base_name, "_", i) }
    roots[[name]] <<- path
  }

  add_root(Sys.getenv("SHINY_ROOT_NAME", "home"),
           Sys.getenv("SHINY_ROOT_PATH", "/browse"))
  add_root("output", "/out")

  if (dir.exists("/mounts")) {
    for (d in list.dirs("/mounts", recursive = FALSE, full.names = TRUE)) {
      add_root(basename(d), d)
    }
  }

  extra <- Sys.getenv("SHINY_EXTRA_ROOTS", "")
  if (nzchar(extra)) {
    for (pair in strsplit(extra, ",", fixed = TRUE)[[1]]) {
      kv <- strsplit(trimws(pair), "=", fixed = TRUE)[[1]]
      if (length(kv) == 2) add_root(trimws(kv[1]), trimws(kv[2]))
    }
  }

  if (length(roots) == 0) roots <- c(root = "/")
  roots
}

server <- function(input, output, session) {
  # roots <- c(home = "~")
  roots <- build_shiny_roots()
  shinyDirChoose(input, "dir", roots = roots)

  mod_degsea_server("degsea", roots = roots)
  mod_signature_proj_server("sign_proj", roots = roots)
  mod_anchored_gsea_server("anchored_gsea", roots = roots)
  mod_outrider_server("outrider", roots = roots)
  mod_univariate_cox_server("univariate_cox", roots = roots)
  mod_tximport_server("tximport", roots = roots)
  mod_signatures_comparison_server("sig_comparison", roots = roots)
  mod_ica_server("ica", roots = roots)

  # mod_gsva_server("gsva", roots = roots)
}

shinyApp(ui, server)