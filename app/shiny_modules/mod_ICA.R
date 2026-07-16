# shiny_modules/mod_ICA.R
# Independent Component Analysis of a bulk transcriptomic matrix.
# Part 1: Most Stable Transcriptomic Dimension (choose the number of components).
# Part 2: run stabilized ICA, save A / S, and visualize each plot type.
# Compute lives in R_analysis/ICA.R -> py/py_ICA.py (via reticulate).

mod_ica_ui <- function(id) {
  ns <- NS(id)

  page_sidebar(
    sidebar = sidebar(
      accordion(
        id = ns("ica_acc"),
        open = FALSE,
        multiple = TRUE,

        accordion_panel(
          "Files *",
          fileInput(ns("bulk_file"), "Expression matrix (genes x samples) : *",
                    accept = c(".csv", ".tsv")),
          fileInput(ns("clinic_file"), "Clinic (optional) :", accept = c(".csv", ".tsv")),
          shinyFiles::shinyDirButton(ns("output_dir"), "Output directory : *", "Upload")
        ),

        accordion_panel(
          "Preprocessing",
          checkboxInput(ns("do_vst"), "VST-normalize (normVST_bulk)", value = TRUE)
        ),

        accordion_panel(
          "1) Most stable dimension",
          numericInput(ns("mstd_min"), "Min number of components :", value = 2, min = 2, step = 1),
          numericInput(ns("mstd_max"), "Max number of components :", value = 40, min = 3, step = 1),
          numericInput(ns("mstd_step"), "Step :", value = 2, min = 1, step = 1),
          numericInput(ns("mstd_runs"), "ICA runs per dimension :", value = 10, min = 2, step = 1),
          actionButton(ns("run_mstd"), "Run MSTD")
        ),

        accordion_panel(
          "2) Run ICA",
          numericInput(ns("n_components"), "Number of components (chosen dim) :",
                       value = 10, min = 2, step = 1),
          numericInput(ns("ica_runs"), "ICA runs (stability) :", value = 50, min = 2, step = 1),
          actionButton(ns("run_ica"), "Run ICA")
        )
      )
    ),

    navset_card_tab(
      nav_panel("Tutorial", ICA_tutorial_card()),
      nav_panel("Logs", card(verbatimTextOutput(ns("logs")))),
      nav_panel("MSTD",
                card(div(class = "responsive-plot-frame",
                         imageOutput(ns("mstd_plot"), width = "100%", height = "100%")))),
      nav_panel("Scatter with marginals",
                card(
                  layout_columns(
                    col_widths = c(3, 3, 3, 3),
                    selectInput(ns("scatter_x"), "Metagene X :", choices = character(0),
                                multiple = FALSE),
                    selectInput(ns("scatter_y"), "Metagene Y :", choices = character(0),
                                multiple = FALSE),
                    selectInput(ns("color_col"), "Color by (clinic) :",
                                choices = c("(none)" = ""), multiple = FALSE),
                    numericInput(ns("scatter_bins"), "Marginal bins :", value = 15, min = 5, step = 1)
                  ),
                  div(class = "responsive-plot-frame",
                      imageOutput(ns("scatter"), width = "100%", height = "100%"))
                )),
      nav_panel("Activity heatmap",
                card(div(class = "responsive-plot-frame",
                         imageOutput(ns("a_heatmap"), width = "100%", height = "100%")))),
      nav_panel("Metagene gene weights",
                card(div(class = "responsive-plot-frame",
                         imageOutput(ns("s_dist"), width = "100%", height = "100%")))),
      nav_panel("Metagene correlation",
                card(div(class = "responsive-plot-frame",
                         imageOutput(ns("corr"), width = "100%", height = "100%")))),
      nav_panel("Clinical association",
                layout_columns(
                  col_width = c(6, 6),
                  card(div(class = "responsive-plot-frame",
                           imageOutput(ns("clinical_continuous"), width = "100%", height = "100%"))),
                  card(div(class = "responsive-plot-frame",
                           imageOutput(ns("clinical_categorical"), width = "100%", height = "100%")))
                )),
      nav_panel("Matrices",
                card(
                  tags$p(class = "text-muted",
                         "Saved under ", tags$code("output_dir/ICA/"),
                         " as CSV + RDS (notebook convention): ",
                         tags$code("A_matrix_*"), " = gene weights (metagenes x genes), ",
                         tags$code("S_matrix_*"), " = sample activities (metagenes x samples)."),
                  tags$strong("A - gene weights (genes x metagenes)"),
                  DT::DTOutput(ns("A_table")),
                  tags$hr(),
                  tags$strong("S - sample activities (samples x metagenes)"),
                  DT::DTOutput(ns("S_table"))
                ))
    )
  )
}


mod_ica_server <- function(id, roots = c(home = "~")) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    library(reticulate)
    Sys.unsetenv("RETICULATE_PYTHON")
    use_python("/opt/conda/envs/BulkTools/bin/python", required = TRUE)
    source_python("py/py_ICA.py")

    bulk_df <- reactive({
      read_delim_auto(input$bulk_file)
    })

    .clinic_in <- clinic_input(reactive(input$clinic_file))
    clinic_df   <- .clinic_in$df
    clinic_note <- .clinic_in$message

    shinyFiles::shinyDirChoose(
      input, id = "output_dir", session = session,
      roots = roots, filetypes = c("", "txt", "tsv", "csv")
    )

    output_dir_path <- reactive({
      req(input$output_dir)
      shinyFiles::parseDirPath(roots, input$output_dir)
    })

    # Metagene pickers follow the chosen number of components (Metagene0..K-1).
    observeEvent(input$n_components, {
      req(input$n_components, input$n_components >= 1)
      comps <- paste0("Metagene", seq_len(as.integer(input$n_components)) - 1)
      updateSelectInput(session, "scatter_x", choices = comps, selected = comps[1])
      updateSelectInput(session, "scatter_y", choices = comps,
                        selected = if (length(comps) >= 2) comps[2] else comps[1])
    })

    # Color-by choices come from the (optional) clinical table.
    observeEvent(input$clinic_file, {
      req(clinic_df())
      updateSelectInput(session, "color_col",
                        choices = c("(none)" = "", names(clinic_df())), selected = "")
    })

    # ---- Part 1: MSTD ----
    mstd_res <- eventReactive(input$run_mstd, {
      req(bulk_df(), output_dir_path())
      withProgress(message = "Estimating stable dimension (MSTD)...", value = 0.05, {
        res <- mstd_pipe(
          rnafilt_counts = bulk_df(),
          output_dir     = path.expand(output_dir_path()),
          m              = input$mstd_min,
          M              = input$mstd_max,
          step           = input$mstd_step,
          n_runs         = input$mstd_runs,
          do_vst         = isTRUE(input$do_vst),
          progress_cb    = function(frac, detail) setProgress(value = frac, detail = detail)
        )
        setProgress(value = 1, detail = "done")
        res
      })
    })

    # ---- Part 2: ICA ----
    ica_res <- eventReactive(input$run_ica, {
      req(bulk_df(), output_dir_path())
      clinic <- if (isTruthy(input$clinic_file)) clinic_df() else NULL
      withProgress(message = "Running stabilized ICA...", value = 0.05, {
        res <- ica_pipe(
          rnafilt_counts = bulk_df(),
          output_dir     = path.expand(output_dir_path()),
          n_components   = input$n_components,
          n_runs         = input$ica_runs,
          do_vst         = isTRUE(input$do_vst),
          clinic_annot   = clinic,
          progress_cb    = function(frac, detail) setProgress(value = frac, detail = detail)
        )
        setProgress(value = 1, detail = "done")
        res
      })
    })

    # Sync the scatter pickers to the actual metagenes once ICA has run,
    # preserving the current selection when still valid.
    observeEvent(ica_res(), {
      comps <- ica_res()$comp_names
      sel_x <- if (isTruthy(input$scatter_x) && input$scatter_x %in% comps) input$scatter_x else comps[1]
      sel_y <- if (isTruthy(input$scatter_y) && input$scatter_y %in% comps) {
        input$scatter_y
      } else if (length(comps) >= 2) comps[2] else comps[1]
      updateSelectInput(session, "scatter_x", choices = comps, selected = sel_x)
      updateSelectInput(session, "scatter_y", choices = comps, selected = sel_y)
    })

    # ---- Scatter with marginals: re-plotted live (no ICA re-run needed) ----
    scatter_path <- reactive({
      req(ica_res(), input$scatter_x, input$scatter_y)
      clinic <- if (isTruthy(input$clinic_file)) clinic_df() else NULL
      plot_scatter_marginals(
        S            = ica_res()$S,
        comp_names   = ica_res()$comp_names,
        sample_names = ica_res()$sample_names,
        comp_x       = input$scatter_x,
        comp_y       = input$scatter_y,
        out_dir      = ica_res()$out_dir,
        clinic       = clinic,
        color_col    = input$color_col,
        bins         = as.integer(input$scatter_bins)
      )
    })

    # ---- Logs ----
    output$logs <- renderPrint({
      msg <- character(0)
      note <- clinic_note()
      if (!is.null(note)) msg <- c(msg, note, "")

      # shiny.silent.error = not run yet / inputs incomplete; a real error means the
      # pipeline failed and must be surfaced here rather than swallowed.
      grab <- function(r) tryCatch(
        r(),
        shiny.silent.error = function(e) NULL,
        error = function(e) e
      )

      m <- grab(mstd_res)
      if (inherits(m, "error")) {
        msg <- c(msg, "MSTD ERROR:", conditionMessage(m), "")
      } else if (!is.null(m)) {
        msg <- c(msg, "MSTD outputs:", paste0("  ", names(m)), "")
      }

      i <- grab(ica_res)
      if (inherits(i, "error")) {
        msg <- c(msg, "ICA ERROR:", conditionMessage(i))
      } else if (!is.null(i)) {
        msg <- c(msg, "ICA outputs:", paste0("  ", names(i)))
      }

      if (length(msg) == 0) msg <- "Run MSTD and/or ICA to populate the tabs."
      cat(msg, sep = "\n")
    })
    # Keep the logs computed even when another tab (e.g. Tutorial) is shown, so that
    # clicking "Run" starts the analysis whatever the active tab is.
    outputOptions(output, "logs", suspendWhenHidden = FALSE)

    # ---- Image outputs ----
    render_png <- function(path_reactive, alt) {
      renderImage({
        p <- path_reactive()
        req(p, file.exists(p))
        list(src = p, contentType = "image/png", alt = alt)
      }, deleteFile = FALSE)
    }

    output$mstd_plot  <- render_png(reactive(mstd_res()$mstd_plot), "MSTD plot")
    output$scatter    <- render_png(scatter_path, "Scatter with marginals")
    output$a_heatmap  <- render_png(reactive(ica_res()$activity_heatmap), "Activity heatmap")
    output$s_dist     <- render_png(reactive(ica_res()$source_dist), "Metagene gene weights")
    output$corr       <- render_png(reactive(ica_res()$corr), "Metagene correlation")
    output$clinical_continuous  <- render_png(reactive(ica_res()$clinical$continuous),
                                              "Clinical association (continuous)")
    output$clinical_categorical <- render_png(reactive(ica_res()$clinical$categorical),
                                              "Clinical association (categorical)")

    # ---- Tables ----
    output$A_table <- DT::renderDT({
      req(ica_res()$A)
      # A is metagenes x genes; transpose so DT paginates over rows (genes).
      A_disp <- as.data.frame(t(as.matrix(ica_res()$A)))
      DT::datatable(round(A_disp, 4), options = list(scrollX = TRUE, pageLength = 10))
    })

    output$S_table <- DT::renderDT({
      req(ica_res()$S)
      # S is metagenes x samples; transpose so DT paginates over rows (samples).
      S_disp <- as.data.frame(t(as.matrix(ica_res()$S)))
      DT::datatable(round(S_disp, 4), options = list(scrollX = TRUE, pageLength = 10))
    })
  })
}
