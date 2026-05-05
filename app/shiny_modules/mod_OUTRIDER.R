# shiny_modules/mod_OUTRIDER.R
mod_outrider_ui <- function(id) {
  ns <- NS(id)

  page_sidebar(
    sidebar = sidebar(
      accordion(
        id = ns("outrider_acc"),
        open = FALSE,
        multiple = TRUE,

        accordion_panel(
          "Files *",
          fileInput(ns("bulk_file"),   "Bulk path : *",   accept = c(".csv", ".tsv")),
          fileInput(ns("clinic_file"), "Clinic path : *", accept = c(".csv", ".tsv")),
          shinyFiles::shinyDirButton(ns("output_dir"), "Output directory : *", "Upload")
        ),

        accordion_panel(
          "Parameters *",
          div(
            selectizeInput(
              ns("samples_to_exclude"), "Samples to exclude :",
              choices = character(0),
              multiple = TRUE,
              options = list(placeholder = "Select sample IDs to exclude…")
            ),
            div(
              class = "d-flex justify-content-end mb-2",
              actionButton(
                ns("clear_samples_to_exclude"),
                label = NULL,
                icon = icon("trash"),
                class = "btn btn-sm btn-outline-danger"
              )
            )
          ),
          div(
            selectInput(
              ns("confounders"), "Confounders :",
              choices = character(0),
              multiple = TRUE, selectize = FALSE, size = 7
            ),
            div(
              class = "d-flex justify-content-end mb-2",
              actionButton(
                ns("clear_confounders"),
                label = NULL,
                icon = icon("trash"),
                class = "btn btn-sm btn-outline-danger"
              )
            )
          ),
          div(
            selectizeInput(
              ns("label_col"), "Grouping column :",
              choices = character(0),
              multiple = FALSE,
              options = list(placeholder = "Select a column for labelling samples in plots…")
            ),
            div(
              class = "d-flex justify-content-end mb-2",
              actionButton(
                ns("clear_label_col"),
                label = NULL,
                icon = icon("trash"),
                class = "btn btn-sm btn-outline-danger"
              )
            )
          ),
          div(
            selectizeInput(
              ns("volcano_samples"), "Samples to volcano :",
              choices = character(0),
              multiple = TRUE,
              options = list(placeholder = "Select sample IDs for volcano plots…")
            ),
            div(
              class = "d-flex justify-content-end mb-2",
              actionButton(
                ns("clear_volcano_samples"),
                label = NULL,
                icon = icon("trash"),
                class = "btn btn-sm btn-outline-danger"
              )
            )
          ),
          div(
            selectizeInput(
              ns("plot_genes"), "Genes to plot :",
              choices = character(0),
              multiple = TRUE,
              options = list(placeholder = "Type a gene symbol…")
            ),
            div(
              class = "d-flex justify-content-end mb-2",
              actionButton(
                ns("clear_plot_genes"),
                label = NULL,
                icon = icon("trash"),
                class = "btn btn-sm btn-outline-danger"
              )
            )
          ),
          numericInput(
            ns("iterations"), "Iterations for confounder control :",
            value = 3, min = 1, max = 20, step = 1
          )
        )
      ),
      actionButton(ns("run_OUTRIDER"), "Run OUTRIDER")
    ),

    navset_card_tab(
      nav_panel("Logs", card(verbatimTextOutput(ns("logs")))),
      nav_panel("Aberrant per sample", card(plotOutput(ns("aberrant_per_sample"), height = 450))),
      nav_panel("Results", card(DT::DTOutput(ns("results_table")))),
      nav_panel(
        "Volcano",
        card(
          div(
            class = "mb-2",
            shinyFiles::shinyFilesButton(
              ns("pick_volcano"),
              label = "Open volcano plot",
              title = "Select a volcano plot to display",
              multiple = FALSE,
              icon = icon("folder-open")
            )
          ),
          tabsetPanel(id = ns("volcano_tabs"))
        )
      ),
      nav_panel(
        "Expression rank",
        card(
          div(
            class = "mb-2",
            shinyFiles::shinyFilesButton(
              ns("pick_rank"),
              label = "Open expression rank plot",
              title = "Select an expression rank plot to display",
              multiple = FALSE,
              icon = icon("folder-open")
            )
          ),
          tabsetPanel(id = ns("rank_tabs"))
        )
      ),
      nav_panel(
        "Expected vs observed",
        card(
          div(
            class = "mb-2",
            shinyFiles::shinyFilesButton(
              ns("pick_exp"),
              label = "Open expected vs observed plot",
              title = "Select an expected vs observed plot to display",
              multiple = FALSE,
              icon = icon("folder-open")
            )
          ),
          tabsetPanel(id = ns("exp_tabs"))
        )
      )
    )
  )
}


mod_outrider_server <- function(id, roots = c(home = "~")) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    clinic_df <- reactive({
      read_delim_auto(input$clinic_file)
    })

    bulk_df <- reactive({
      read_delim_auto(input$bulk_file)
    })

    sample_id_choices <- reactive({
      ids <- character(0)
      if (isTruthy(input$bulk_file)) {
        ids <- colnames(bulk_df())
      }
      if (isTruthy(input$clinic_file)) {
        clinic <- clinic_df()
        clinic_ids <- if ("ID_Patient" %in% colnames(clinic)) {
          as.character(clinic$ID_Patient)
        } else {
          rownames(clinic)
        }
        ids <- if (length(ids) == 0) clinic_ids else intersect(ids, clinic_ids)
      }
      sort(unique(ids[!is.na(ids) & nzchar(ids)]))
    })

    observeEvent(input$bulk_file, {
      req(bulk_df())
      genes <- rownames(bulk_df())
      updateSelectizeInput(session, "plot_genes", choices = genes, server = TRUE)
      updateSelectizeInput(session, "samples_to_exclude", choices = sample_id_choices(), server = TRUE)
      updateSelectizeInput(session, "volcano_samples",   choices = sample_id_choices(), server = TRUE)
    })

    observeEvent(input$clinic_file, {
      req(clinic_df())
      updateSelectInput(session, "label_col", choices = names(clinic_df()), selected = NULL)
      updateSelectInput(session, "confounders", choices = names(clinic_df()), selected = NULL)
      updateSelectizeInput(session, "samples_to_exclude", choices = sample_id_choices(), server = TRUE)
      updateSelectizeInput(session, "volcano_samples",   choices = sample_id_choices(), server = TRUE)
    })


    shinyFiles::shinyDirChoose(
      input,
      id = "output_dir",
      session = session,
      roots = roots,
      filetypes = c("", "txt", "bigWig", "tsv", "csv", "bw")
    )

    output_dir_path <- reactive({
      req(input$output_dir)
      shinyFiles::parseDirPath(roots, input$output_dir)
    })

    # ---- Per-plot browsers (volcano / expression rank / expected vs observed) ----
    plot_dir_roots <- function(plot_subdir) {
      reactive({
        out <- tryCatch(output_dir_path(), error = function(e) NULL)
        if (is.null(out) || !nzchar(out)) {
          return(roots)
        }
        plot_dir <- file.path(out, "OUTRIDER", "plots", plot_subdir)
        starting <- if (dir.exists(plot_dir)) plot_dir else out
        c(plots = starting, output = out, home = roots[["home"]])
      })
    }

    setup_plot_picker <- function(file_input_id, tabs_id, plot_subdir, output_prefix) {
      roots_reactive <- plot_dir_roots(plot_subdir)
      state <- reactiveValues(open_paths = character(0), counter = 0L)

      shinyFiles::shinyFileChoose(
        input,
        id = file_input_id,
        session = session,
        roots = roots_reactive,
        filetypes = c("png", "jpg", "jpeg", "svg")
      )

      observeEvent(input[[file_input_id]], {
        sel <- input[[file_input_id]]
        if (is.null(sel) || is.integer(sel)) return()

        parsed <- shinyFiles::parseFilePaths(roots_reactive(), sel)
        if (nrow(parsed) == 0) return()

        file_path <- as.character(parsed$datapath[1])

        if (file_path %in% state$open_paths) {
          updateTabsetPanel(session, tabs_id, selected = file_path)
          return()
        }

        state$counter <- state$counter + 1L
        output_id <- paste0(output_prefix, "_img_", state$counter)

        local({
          img_path <- file_path
          ext <- tolower(tools::file_ext(img_path))
          content_type <- switch(
            ext,
            png = "image/png",
            jpg = "image/jpeg",
            jpeg = "image/jpeg",
            svg = "image/svg+xml",
            "image/png"
          )
          output[[output_id]] <- renderImage({
            list(src = img_path, contentType = content_type, alt = basename(img_path))
          }, deleteFile = FALSE)
        })

        insertTab(
          inputId = tabs_id,
          tabPanel(
            title = basename(file_path),
            value = file_path,
            imageOutput(ns(output_id), height = "600px")
          ),
          session = session,
          select = TRUE
        )

        state$open_paths <- c(state$open_paths, file_path)
      })
    }

    setup_plot_picker("pick_volcano", "volcano_tabs", "volcanoes",            "volcano")
    setup_plot_picker("pick_rank",    "rank_tabs",    "expression_rank",      "rank")
    setup_plot_picker("pick_exp",     "exp_tabs",     "expected_vs_observed", "exp")

    # Clear buttons for selection fields
    observeEvent(input$clear_samples_to_exclude, {
      updateSelectizeInput(session, "samples_to_exclude", selected = character(0))
    })
    observeEvent(input$clear_confounders, {
      updateSelectInput(session, "confounders", selected = character(0))
    })
    observeEvent(input$clear_volcano_samples, {
      updateSelectizeInput(session, "volcano_samples", selected = character(0))
    })
    observeEvent(input$clear_plot_genes, {
      updateSelectizeInput(session, "plot_genes", selected = character(0))
    })

    # Empty multi-selects normalize to NULL before being passed to the pipe
    or_NULL <- function(x) {
      if (is.null(x) || length(x) == 0 || all(!nzchar(x))) NULL else x
    }

    # ---- RUN OUTRIDER (this must NOT be inside observeEvent) ----
    outrider_res <- eventReactive(input$run_OUTRIDER, {
      req(
        clinic_df(), bulk_df(),
        output_dir_path(), input$iterations
      )

      withProgress(message = "Running OUTRIDER...", value = 0, {
        incProgress(0.1)
        res <- OUTRIDER_pipe(
          rnafilt_counts     = bulk_df(),
          clinic_annot       = clinic_df(),
          output_dir         = output_dir_path(),
          samples_to_exclude = or_NULL(input$samples_to_exclude),
          confounders        = or_NULL(input$confounders),
          volcano_samples    = or_NULL(input$volcano_samples),
          plot_genes         = or_NULL(input$plot_genes),
          label_column       = or_NULL(input$label_col),
          iterations         = input$iterations
        )
        incProgress(1)
        res
      })
    })

    # outputs
    output$logs <- renderPrint({
      req(outrider_res())
      names(outrider_res())
    })

    output$aberrant_per_sample <- renderPlot({
      p <- outrider_res()$plots$aberrant_per_sample
      req(p)
      print(p)
    })

    output$results_table <- DT::renderDT({
      df <- outrider_res()$tables$results
      req(df)
      df
    })
  })
}
