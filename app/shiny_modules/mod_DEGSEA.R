# shiny_modules/mod_DEGSEA.R
mod_degsea_ui <- function(id) {
  ns <- NS(id)

  page_sidebar(
    sidebar = sidebar(
      accordion(
        id = ns("degsea_acc"),
        open = FALSE,
        multiple = TRUE,

        accordion_panel(
          "Files *",
          fileInput(ns("bulk_file"),   "Bulk path : *",   accept = c(".csv", ".tsv")),
          fileInput(ns("clinic_file"), "Clinic path : *", accept = c(".csv", ".tsv")),
          shinyFiles::shinyDirButton(ns("output_dir"), "Output directory : *", "Upload")
        ),

        accordion_panel(
          "Filtering on clinic",
          div(
            class = "d-flex justify-content-between align-items-center mb-2",
            actionButton(
              ns("add_clinic_filter"),
              label = "Add clinical filter",
              icon = icon("plus"),
              class = "btn btn-sm btn-outline-primary"
            )
          ),
          uiOutput(ns("clinic_filters_ui"))
        ),

        accordion_panel(
          "Filtering on gene",
          selectizeInput(
            ns("gene_filt"), "Filter on gene :",
            choices = stats::setNames("", ""),
            selected = "",
            multiple = FALSE,
            options = list(placeholder = "Type a gene symbol…")
          ),
          numericInput(
            ns("quantile"), "Quantile threshold :",
            value = 0, min = 0, max = 1, step = 0.05
          ),
          selectInput(
            ns("low_or_high"), "Keep low or high expression :",
            choices = c("low", "high"),
            multiple = FALSE, selectize = FALSE, size = 2
          )
        ),

        accordion_panel(
          "Parameters *",
          selectInput(
            ns("DESeq_covar"), "DESeq covariates :",
            choices = character(0),
            multiple = TRUE, selectize = FALSE, size = 7
          ),
          div(
            class = "border rounded p-3 mb-3",
            div(
              class = "d-flex justify-content-between align-items-center mb-2",
              actionButton(
                ns("add_control_filter"),
                label = "Select CONTROL group on clinical columns",
                icon = icon("plus"),
                class = "btn btn-sm btn-outline-primary"
              )
            ),
            uiOutput(ns("control_filters_ui")),
            div(
              class = "d-flex justify-content-between align-items-center mb-2",
              actionButton(
                ns("add_control_gene_filter"),
                label = "Select CONTROL group on gene expression",
                icon = icon("plus"),
                class = "btn btn-sm btn-outline-primary"
              )
            ),
            uiOutput(ns("control_gene_filters_ui"))
          ),
          div(
            class = "border rounded p-3 mb-3",
            div(
              class = "d-flex justify-content-between align-items-center mb-2",
              actionButton(
                ns("add_test_filter"),
                label = "Select TEST group on clinical columns",
                icon = icon("plus"),
                class = "btn btn-sm btn-outline-primary"
              )
            ),
            uiOutput(ns("test_filters_ui")),
            div(
              class = "d-flex justify-content-between align-items-center mb-2",
              actionButton(
                ns("add_test_gene_filter"),
                label = "Select TEST group on gene expression",
                icon = icon("plus"),
                class = "btn btn-sm btn-outline-primary"
              )
            ),
            uiOutput(ns("test_gene_filters_ui"))
          ),
          selectInput(
            ns("GSEA_geneset"), "GSEA geneset : *",
            choices = c("all_pathways", "REACTOME_pathways", "GOBP_pathways", 
                        "KEGG_pathways", "QoL_pathways", "transcript_factor_target",
                        "boyault_sets"),
            multiple = FALSE, selectize = FALSE, size = 4
          )
        )
      ),
      actionButton(ns("run_DEGSEA"), "Run DEGSEA")
    ),

    navset_card_tab(
      nav_panel("Tutorial", degsea_tutorial_card()),
      nav_panel("Logs", card(verbatimTextOutput(ns("logs")))),
      nav_panel("Volcano", card(plotOutput(ns("volcano"), height = 450))),
      nav_panel("Differential Expression", card(DT::DTOutput(ns("top_DE")))),
      nav_panel("GSEA table",  card(DT::DTOutput(ns("top_GSEA")))),
      nav_panel("GESECA table", card(DT::DTOutput(ns("top_GESECA")))),
      nav_panel("ssGSEA table", card(DT::DTOutput(ns("top_ssGSEA"))))
    )
  )
}

mod_degsea_server <- function(id, roots = c(home = "~")) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    clinic_df <- reactive({
      read_delim_auto(input$clinic_file)
    })

    bulk_df <- reactive({
      read_delim_auto(input$bulk_file)
    })

    gene_filter_choices <- function(genes = character(0)) {
      c(stats::setNames("", ""), stats::setNames(genes, genes))
    }

    clinic_modality_choices <- function(clinic_col) {
      req(clinic_df(), clinic_col)
      values <- clinic_df()[[clinic_col]]
      choices <- if (is.factor(values)) levels(values) else sort(unique(as.character(values)))
      choices[!is.na(choices) & nzchar(choices)]
    }

    create_filter_state <- function() {
      reactiveValues(ids = integer(), next_id = 0)
    }

    reset_filter_state <- function(filter_state) {
      filter_state$ids <- integer()
    }

    collect_named_filters <- function(filter_state, filter_prefix, filter_label) {
      reactive({
        if (length(filter_state$ids) == 0) {
          return(NULL)
        }

        selected_filters <- list()

        for (filter_id in filter_state$ids) {
          clinic_col <- input[[paste0(filter_prefix, "_filter_col_", filter_id)]]
          clinic_col <- trimws(if (is.null(clinic_col)) "" else as.character(clinic_col))

          modalities <- input[[paste0(filter_prefix, "_filter_modalities_", filter_id)]]
          modalities <- as.character(modalities)
          modalities <- trimws(modalities)
          modalities <- modalities[!is.na(modalities) & nzchar(modalities)]

          if (!nzchar(clinic_col)) {
            next
          }

          if (length(modalities) == 0) {
            stop(sprintf("Please select at least one modality for the %s filter '%s'.", filter_label, clinic_col))
          }

          selected_filters[[clinic_col]] <- sort(unique(c(selected_filters[[clinic_col]], modalities)))
        }

        if (length(selected_filters) == 0) {
          return(NULL)
        }

        selected_filters[sort(names(selected_filters))]
      })
    }

    collect_group_gene_filters <- function(filter_state, filter_prefix) {
      reactive({
        if (length(filter_state$ids) == 0) {
          return(NULL)
        }

        selected_filters <- list()

        for (filter_id in filter_state$ids) {
          gene_name <- input[[paste0(filter_prefix, "_gene_filter_gene_", filter_id)]]
          gene_name <- trimws(if (is.null(gene_name)) "" else as.character(gene_name))

          if (!nzchar(gene_name)) {
            next
          }

          selected_filters[[length(selected_filters) + 1]] <- list(
            gene = gene_name,
            quantile_thr = input[[paste0(filter_prefix, "_gene_filter_quantile_", filter_id)]],
            keep_low_or_high = input[[paste0(filter_prefix, "_gene_filter_direction_", filter_id)]]
          )
        }

        normalize_group_gene_filters(selected_filters)
      })
    }

    render_filter_set <- function(filter_state, filter_prefix, empty_message=NULL) {
      output[[paste0(filter_prefix, "_filters_ui")]] <- renderUI({
        clinic_cols <- if (isTruthy(input$clinic_file)) names(clinic_df()) else character(0)

        if (length(filter_state$ids) == 0) {
          return(tags$p(class = "text-muted mb-0", empty_message))
        }

        tagList(lapply(filter_state$ids, function(filter_id) {
          col_input_id <- paste0(filter_prefix, "_filter_col_", filter_id)
          mod_input_id <- paste0(filter_prefix, "_filter_modalities_", filter_id)
          remove_input_id <- paste0("remove_", filter_prefix, "_filter_", filter_id)

          selected_col <- isolate(input[[col_input_id]])
          selected_modalities <- isolate(input[[mod_input_id]])
          modality_choices <- if (!is.null(selected_col) && nzchar(selected_col)) {
            clinic_modality_choices(selected_col)
          } else {
            character(0)
          }
          selected_modalities <- selected_modalities[selected_modalities %in% modality_choices]

          div(
            class = "border rounded p-3 mb-2",
            selectInput(
              ns(col_input_id), "Filter on clin column :",
              choices = clinic_cols,
              selected = selected_col,
              multiple = FALSE, selectize = FALSE, size = 7
            ),
            selectInput(
              ns(mod_input_id), "Kept modality :",
              choices = modality_choices,
              selected = selected_modalities,
              multiple = TRUE, selectize = FALSE, size = 7
            ),
            div(
              class = "d-flex justify-content-end",
              actionButton(
                ns(remove_input_id),
                label = NULL,
                icon = icon("trash"),
                class = "btn btn-sm btn-outline-danger"
              )
            )
          )
        }))
      })
    }

    render_gene_filter_set <- function(filter_state, filter_prefix, empty_message=NULL) {
      output[[paste0(filter_prefix, "_gene_filters_ui")]] <- renderUI({
        gene_choices <- if (isTruthy(input$bulk_file)) {
          gene_filter_choices(rownames(bulk_df()))
        } else {
          gene_filter_choices()
        }

        if (length(filter_state$ids) == 0) {
          return(tags$p(class = "text-muted mb-0", empty_message))
        }

        tagList(lapply(filter_state$ids, function(filter_id) {
          gene_input_id <- paste0(filter_prefix, "_gene_filter_gene_", filter_id)
          quantile_input_id <- paste0(filter_prefix, "_gene_filter_quantile_", filter_id)
          direction_input_id <- paste0(filter_prefix, "_gene_filter_direction_", filter_id)
          remove_input_id <- paste0("remove_", filter_prefix, "_gene_filter_", filter_id)

          selected_gene <- isolate(input[[gene_input_id]])
          if (is.null(selected_gene) || !nzchar(selected_gene) || !selected_gene %in% names(gene_choices)) {
            selected_gene <- ""
          }

          selected_quantile <- isolate(input[[quantile_input_id]])
          if (is.null(selected_quantile) || length(selected_quantile) == 0 || is.na(selected_quantile)) {
            selected_quantile <- 0.5
          }

          selected_direction <- isolate(input[[direction_input_id]])
          if (is.null(selected_direction) || !selected_direction %in% c("low", "high")) {
            selected_direction <- "low"
          }

          div(
            class = "border rounded p-3 mb-2",
            selectizeInput(
              ns(gene_input_id), "Expression gene :",
              choices = gene_choices,
              selected = selected_gene,
              multiple = FALSE,
              options = list(placeholder = "Type a gene symbol…")
            ),
            numericInput(
              ns(quantile_input_id), "Expression quantile :",
              value = selected_quantile, min = 0, max = 1, step = 0.05
            ),
            selectInput(
              ns(direction_input_id), "Keep low or high expression :",
              choices = c("low (< quantile)" = "low", "high (> quantile)" = "high"),
              selected = selected_direction,
              multiple = FALSE, selectize = FALSE, size = 2
            ),
            div(
              class = "d-flex justify-content-end",
              actionButton(
                ns(remove_input_id),
                label = NULL,
                icon = icon("trash"),
                class = "btn btn-sm btn-outline-danger"
              )
            )
          )
        }))
      })
    }

    register_filter_set <- function(filter_state, filter_prefix) {
      observeEvent(input[[paste0("add_", filter_prefix, "_filter")]], {
        filter_id <- filter_state$next_id + 1
        filter_state$next_id <- filter_id
        filter_state$ids <- c(filter_state$ids, filter_id)

        local({
          current_filter_id <- filter_id
          current_col_id <- paste0(filter_prefix, "_filter_col_", current_filter_id)
          current_mod_id <- paste0(filter_prefix, "_filter_modalities_", current_filter_id)
          current_remove_id <- paste0("remove_", filter_prefix, "_filter_", current_filter_id)

          observeEvent(input[[current_col_id]], {
            selected_col <- input[[current_col_id]]
            modality_choices <- if (!is.null(selected_col) && nzchar(selected_col)) {
              clinic_modality_choices(selected_col)
            } else {
              character(0)
            }

            selected_modalities <- input[[current_mod_id]]
            selected_modalities <- selected_modalities[selected_modalities %in% modality_choices]

            updateSelectInput(
              session,
              current_mod_id,
              choices = modality_choices,
              selected = selected_modalities
            )
          }, ignoreInit = TRUE)

          observeEvent(input[[current_remove_id]], {
            filter_state$ids <- setdiff(filter_state$ids, current_filter_id)
          }, ignoreInit = TRUE)
        })
      })
    }

    register_gene_filter_set <- function(filter_state, filter_prefix) {
      observeEvent(input[[paste0("add_", filter_prefix, "_gene_filter")]], {
        filter_id <- filter_state$next_id + 1
        filter_state$next_id <- filter_id
        filter_state$ids <- c(filter_state$ids, filter_id)

        local({
          current_filter_id <- filter_id
          current_remove_id <- paste0("remove_", filter_prefix, "_gene_filter_", current_filter_id)

          observeEvent(input[[current_remove_id]], {
            filter_state$ids <- setdiff(filter_state$ids, current_filter_id)
          }, ignoreInit = TRUE)
        })
      })
    }

    clinic_filter_state <- create_filter_state()
    control_filter_state <- create_filter_state()
    test_filter_state <- create_filter_state()
    control_gene_filter_state <- create_filter_state()
    test_gene_filter_state <- create_filter_state()

    clinic_filters <- collect_named_filters(clinic_filter_state, "clinic", "clinical")
    control_filters <- collect_named_filters(control_filter_state, "control", "control group")
    test_filters <- collect_named_filters(test_filter_state, "test", "test group")
    control_gene_filter <- collect_group_gene_filters(control_gene_filter_state, "control")
    test_gene_filter <- collect_group_gene_filters(test_gene_filter_state, "test")

    render_filter_set(clinic_filter_state, "clinic")
    render_filter_set(control_filter_state, "control")
    render_filter_set(test_filter_state, "test")
    render_gene_filter_set(control_gene_filter_state, "control")
    render_gene_filter_set(test_gene_filter_state, "test")

    register_filter_set(clinic_filter_state, "clinic")
    register_filter_set(control_filter_state, "control")
    register_filter_set(test_filter_state, "test")
    register_gene_filter_set(control_gene_filter_state, "control")
    register_gene_filter_set(test_gene_filter_state, "test")

    # update gene list when bulk changes
    observeEvent(input$bulk_file, {
      req(bulk_df())
      genes <- rownames(bulk_df())
      gene_choices <- gene_filter_choices(genes)
      updateSelectizeInput(session, "gene_filt", choices = gene_choices, selected = "", server = TRUE)
      reset_filter_state(control_gene_filter_state)
      reset_filter_state(test_gene_filter_state)
    })

    # update clinic column lists
    observeEvent(input$clinic_file, {
      req(clinic_df())
      updateSelectInput(session, "DESeq_covar",    choices = names(clinic_df()), selected = NULL)
      reset_filter_state(clinic_filter_state)
      reset_filter_state(control_filter_state)
      reset_filter_state(test_filter_state)
    })

    # shinyFiles directory chooser (IMPORTANT: pass session=)
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

    # ---- RUN DEGSEA (this must NOT be inside observeEvent) ----
    degsea_res <- eventReactive(input$run_DEGSEA, {
      req(
        clinic_df(), bulk_df(),
        input$GSEA_geneset, output_dir_path()
      )

      withProgress(message = "Running DEGSEA...", value = 0, {
        incProgress(0.2)
        res <- DEGSEA_pipe(
          rnafilt_counts     = bulk_df(),
          clinic_annot       = clinic_df(),
          control_filters    = control_filters(),
          test_filters       = test_filters(),
          pathways_to_use    = input$GSEA_geneset,
          output_dir         = output_dir_path(),
          covariates         = input$DESeq_covar,
          clinic_filters     = clinic_filters(),
          control_gene_filter = control_gene_filter(),
          test_gene_filter    = test_gene_filter(),
          filter_by_gene     = input$gene_filt,
          quantile_thr       = input$quantile,
          keep_low_or_high   = input$low_or_high
        )
        incProgress(1)
        res
      })
    })

    # outputs
    output$logs <- renderPrint({
      req(degsea_res())
      names(degsea_res())
    })

    output$volcano <- renderPlot({
      p <- degsea_res()$plots$volcano
      req(p)
      print(p)
    })

    output$top_DE <- DT::renderDT({
      df <- degsea_res()$tables$top_DE
      req(df)
      df
    })

    output$top_GSEA <- DT::renderDT({
      df <- degsea_res()$tables$top_GSEA
      req(df)
      df
    })

    output$top_GESECA <- DT::renderDT({
      df <- degsea_res()$tables$top_GESECA
      req(df)
      df
    })

    output$top_ssGSEA <- DT::renderDT({
      df <- degsea_res()$tables$top_ssGSEA
      req(df)
      df
    })
  })
}
