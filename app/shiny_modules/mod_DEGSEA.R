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
            tags$span("Add one or more clinical filters"),
            actionButton(
              ns("add_clinic_filter"),
              label = NULL,
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
            choices = character(0),
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
          selectInput(
            ns("DESeq_contrast"), "DESeq contrast : *",
            choices = character(0),
            multiple = FALSE, selectize = FALSE, size = 7
          ),
          selectInput(
            ns("DESeq_control"), "DESeq control : *",
            choices = character(0),
            multiple = TRUE, selectize = FALSE, size = 5
          ),
          selectInput(
            ns("DESeq_test"), "DESeq test : *",
            choices = character(0),
            multiple = TRUE, selectize = FALSE, size = 5
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

    clinic_filter_state <- reactiveValues(ids = integer(), next_id = 0)

    clinic_modality_choices <- function(clinic_col) {
      req(clinic_df(), clinic_col)
      values <- clinic_df()[[clinic_col]]
      choices <- if (is.factor(values)) levels(values) else sort(unique(as.character(values)))
      choices[!is.na(choices) & nzchar(choices)]
    }

    clinic_filters <- reactive({
      if (length(clinic_filter_state$ids) == 0) {
        return(NULL)
      }

      selected_filters <- list()

      for (filter_id in clinic_filter_state$ids) {
        clinic_col <- input[[paste0("clinic_filter_col_", filter_id)]]
        clinic_col <- trimws(if (is.null(clinic_col)) "" else as.character(clinic_col))

        modalities <- input[[paste0("clinic_filter_modalities_", filter_id)]]
        modalities <- as.character(modalities)
        modalities <- trimws(modalities)
        modalities <- modalities[!is.na(modalities) & nzchar(modalities)]

        if (!nzchar(clinic_col)) {
          next
        }

        if (length(modalities) == 0) {
          stop(sprintf("Please select at least one modality for clinical filter '%s'.", clinic_col))
        }

        selected_filters[[clinic_col]] <- sort(unique(c(selected_filters[[clinic_col]], modalities)))
      }

      if (length(selected_filters) == 0) {
        return(NULL)
      }

      selected_filters[sort(names(selected_filters))]
    })

    output$clinic_filters_ui <- renderUI({
      clinic_cols <- if (isTruthy(input$clinic_file)) names(clinic_df()) else character(0)

      if (length(clinic_filter_state$ids) == 0) {
        return(tags$p(class = "text-muted mb-0", "No clinical filter added yet."))
      }

      tagList(lapply(clinic_filter_state$ids, function(filter_id) {
        col_input_id <- paste0("clinic_filter_col_", filter_id)
        mod_input_id <- paste0("clinic_filter_modalities_", filter_id)
        remove_input_id <- paste0("remove_clinic_filter_", filter_id)

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

    # update gene list when bulk changes
    observeEvent(input$bulk_file, {
      req(bulk_df())
      genes <- rownames(bulk_df())
      updateSelectizeInput(session, "gene_filt", choices = genes, server = TRUE)
    })

    # update clinic column lists
    observeEvent(input$clinic_file, {
      req(clinic_df())
      updateSelectInput(session, "DESeq_contrast", choices = names(clinic_df()), selected = "")
      updateSelectInput(session, "DESeq_covar",    choices = names(clinic_df()), selected = NULL)
      clinic_filter_state$ids <- integer()
    })

    observeEvent(input$add_clinic_filter, {
      filter_id <- clinic_filter_state$next_id + 1
      clinic_filter_state$next_id <- filter_id
      clinic_filter_state$ids <- c(clinic_filter_state$ids, filter_id)

      local({
        current_filter_id <- filter_id
        current_col_id <- paste0("clinic_filter_col_", current_filter_id)
        current_mod_id <- paste0("clinic_filter_modalities_", current_filter_id)
        current_remove_id <- paste0("remove_clinic_filter_", current_filter_id)

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
          clinic_filter_state$ids <- setdiff(clinic_filter_state$ids, current_filter_id)
        }, ignoreInit = TRUE)
      })
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

    # contrast -> control/test modalities
    observeEvent(input$DESeq_contrast, {
      req(clinic_df(), input$DESeq_contrast)
      v <- clinic_df()[[input$DESeq_contrast]]
      mods <- if (is.factor(v)) levels(v) else sort(unique(as.character(v)))
      mods <- mods[!is.na(mods) & mods != ""]
      updateSelectInput(session, "DESeq_control", choices = mods, selected = "")
      updateSelectInput(session, "DESeq_test",    choices = mods, selected = "")
    })

    # ---- RUN DEGSEA (this must NOT be inside observeEvent) ----
    degsea_res <- eventReactive(input$run_DEGSEA, {
      req(
        clinic_df(), bulk_df(),
        input$DESeq_contrast, input$DESeq_control, input$DESeq_test,
        input$GSEA_geneset, output_dir_path()
      )

      withProgress(message = "Running DEGSEA...", value = 0, {
        incProgress(0.2)
        res <- DEGSEA_pipe(
          rnafilt_counts     = bulk_df(),
          clinic_annot       = clinic_df(),
          contrast           = input$DESeq_contrast,
          control            = input$DESeq_control,
          test               = input$DESeq_test,
          pathways_to_use    = input$GSEA_geneset,
          output_dir         = output_dir_path(),
          covariates         = input$DESeq_covar,
          clinic_filters     = clinic_filters(),
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
