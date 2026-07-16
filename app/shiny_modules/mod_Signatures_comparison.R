# shiny_modules/mod_Signatures_comparison.R
library(jsonlite)

mod_signatures_comparison_ui <- function(id) {
  ns <- NS(id)

  page_sidebar(
    sidebar = sidebar(
      accordion(
        id = ns("sigcomp_acc"),
        open = FALSE,
        multiple = TRUE,

        accordion_panel(
          "Files *",
          fileInput(ns("bulk_file"),   "Bulk path : *",   accept = c(".csv", ".tsv")),
          fileInput(ns("clinic_file"), "Clinic path : *", accept = c(".csv", ".tsv")),
          shinyFiles::shinyDirButton(ns("output_dir"), "Output directory : *", "Upload")
        ),

        accordion_panel(
          "Filter on clinic & define conditions",
          div(
            class = "d-flex justify-content-between align-items-center mb-2",
            tags$span("Global clinical filters (applied first)"),
            actionButton(ns("add_clinic_filter"), label = NULL, icon = icon("plus"),
                         class = "btn btn-sm btn-outline-primary")
          ),
          uiOutput(ns("clinic_filters_ui")),

          tags$hr(),
          div(
            class = "d-flex justify-content-between align-items-center mb-2",
            tags$strong("Conditions"),
            actionButton(ns("add_condition"), label = NULL, icon = icon("plus"),
                         class = "btn btn-sm btn-outline-primary")
          ),
          uiOutput(ns("conditions_ui")),

          tags$hr(),
          tags$strong("Response"),
          selectInput(ns("response_col"), "Response column :", choices = character(0),
                      multiple = FALSE, selectize = FALSE, size = 5),
          selectInput(ns("responders"), "Responder modalities :", choices = character(0),
                      multiple = TRUE, selectize = FALSE, size = 5),
          selectInput(ns("non_responders"), "Non-responder modalities :", choices = character(0),
                      multiple = TRUE, selectize = FALSE, size = 5)
        ),

        accordion_panel(
          "Select signatures *",
          div(
            class = "d-flex justify-content-between align-items-center mb-2",
            tags$span("Add signatures to compare"),
            actionButton(ns("add_signature"), label = NULL, icon = icon("plus"),
                         class = "btn btn-sm btn-outline-primary")
          ),
          uiOutput(ns("signatures_ui")),
          tags$hr(),
          tags$strong("Add individual genes as signatures"),
          selectizeInput(ns("genes"), "Genes :",
                         choices = character(0), selected = character(0),
                         multiple = TRUE,
                         options = list(placeholder = "Type gene symbols…")),
          tags$hr(),
          tags$strong("Continuous clinical variables (optional)"),
          selectizeInput(ns("clin_continuous"), "Variables :",
                         choices = character(0), selected = character(0),
                         multiple = TRUE,
                         options = list(placeholder = "Numeric clinic columns…"))
        ),

        accordion_panel(
          "Parameters",
          checkboxInput(ns("do_vst"), "VST-normalize (normVST_bulk)", value = TRUE),
          selectInput(ns("corr_method"), "Correlation test :",
                      choices = c("spearman", "pearson", "kendall"),
                      selected = "spearman", multiple = FALSE, selectize = FALSE, size = 3),
          checkboxInput(ns("corr_fdr"), "Use BH-corrected q-values for stars", value = FALSE),
          textInput(ns("responder_label"), "Responder label :", value = "R"),
          textInput(ns("nonresponder_label"), "Non-responder label :", value = "NR")
        )
      ),
      actionButton(ns("run_sigcomp"), "Run comparison")
    ),

    navset_card_tab(
      nav_panel("Tutorial", signatures_comparison_tutorial_card()),
      nav_panel("Logs", card(verbatimTextOutput(ns("logs")))),
      nav_panel("Boxplots",
                card(div(class = "responsive-plot-frame",
                         imageOutput(ns("boxplot"), width = "100%", height = "100%")))),
      nav_panel("Clinical boxplots",
                card(div(class = "responsive-plot-frame",
                         imageOutput(ns("clinical_boxplot"), width = "100%", height = "100%")))),
      nav_panel("Correlations", uiOutput(ns("corr_tabs_ui"))),
      nav_panel("ROC curves",
                card(div(class = "responsive-plot-frame",
                         imageOutput(ns("roc"), width = "100%", height = "100%")))),
      nav_panel("Wilcoxon stats", card(DT::DTOutput(ns("stats_table"))))
    )
  )
}


mod_signatures_comparison_server <- function(id, roots = c(home = "~")) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    `%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

    .clinic_in <- clinic_input(reactive(input$clinic_file))
    clinic_df   <- .clinic_in$df
    clinic_note <- .clinic_in$message
    bulk_df   <- reactive({ read_delim_auto(input$bulk_file) })

    clinic_modality_choices <- function(clinic_col) {
      req(clinic_df(), clinic_col)
      values <- clinic_df()[[clinic_col]]
      choices <- if (is.factor(values)) levels(values) else sort(unique(as.character(values)))
      choices[!is.na(choices) & nzchar(choices)]
    }

    ## ---- Global clinic filters (reuse the projection-module pattern) ----
    clinic_filter_state <- reactiveValues(ids = integer(), next_id = 0)

    clinic_filters <- reactive({
      if (length(clinic_filter_state$ids) == 0) return(NULL)
      selected_filters <- list()
      for (filter_id in clinic_filter_state$ids) {
        clinic_col <- input[[paste0("clinic_filter_col_", filter_id)]]
        clinic_col <- trimws(if (is.null(clinic_col)) "" else as.character(clinic_col))
        modalities <- input[[paste0("clinic_filter_modalities_", filter_id)]]
        modalities <- trimws(as.character(modalities))
        modalities <- modalities[!is.na(modalities) & nzchar(modalities)]
        if (!nzchar(clinic_col)) next
        if (length(modalities) == 0) {
          stop(sprintf("Please select at least one modality for clinical filter '%s'.", clinic_col))
        }
        selected_filters[[clinic_col]] <- sort(unique(c(selected_filters[[clinic_col]], modalities)))
      }
      if (length(selected_filters) == 0) return(NULL)
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
        } else character(0)
        selected_modalities <- selected_modalities[selected_modalities %in% modality_choices]
        div(
          class = "border rounded p-3 mb-2",
          selectInput(ns(col_input_id), "Filter on clin column :", choices = clinic_cols,
                      selected = selected_col, multiple = FALSE, selectize = FALSE, size = 5),
          selectInput(ns(mod_input_id), "Kept modality :", choices = modality_choices,
                      selected = selected_modalities, multiple = TRUE, selectize = FALSE, size = 5),
          div(class = "d-flex justify-content-end",
              actionButton(ns(remove_input_id), label = NULL, icon = icon("trash"),
                           class = "btn btn-sm btn-outline-danger"))
        )
      }))
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
          } else character(0)
          selected_modalities <- input[[current_mod_id]]
          selected_modalities <- selected_modalities[selected_modalities %in% modality_choices]
          updateSelectInput(session, current_mod_id, choices = modality_choices,
                            selected = selected_modalities)
        }, ignoreInit = TRUE)
        observeEvent(input[[current_remove_id]], {
          clinic_filter_state$ids <- setdiff(clinic_filter_state$ids, current_filter_id)
        }, ignoreInit = TRUE)
      })
    })

    ## ---- Populate response column + continuous-variable choices on clinic load ----
    numeric_clinic_cols <- function(df) {
      names(df)[vapply(df, function(x) {
        if (is.numeric(x)) return(TRUE)
        xn <- suppressWarnings(as.numeric(as.character(x)))
        n_ok  <- sum(!is.na(xn))
        n_val <- sum(!is.na(x) & nzchar(as.character(x)))
        n_val > 0 && n_ok >= 0.8 * n_val && length(unique(xn[!is.na(xn)])) > 2
      }, logical(1))]
    }

    observeEvent(input$clinic_file, {
      req(clinic_df())
      updateSelectInput(session, "response_col", choices = names(clinic_df()), selected = "")
      updateSelectizeInput(session, "clin_continuous",
                           choices = numeric_clinic_cols(clinic_df()),
                           selected = character(0), server = TRUE)
      clinic_filter_state$ids <- integer()
    })

    ## ---- Conditions: dynamic list (>= 1), "+" to add ----
    cond_state <- reactiveValues(ids = integer(), next_id = 0)

    output$conditions_ui <- renderUI({
      clinic_cols <- if (isTruthy(input$clinic_file)) names(clinic_df()) else character(0)
      if (length(cond_state$ids) == 0) {
        return(tags$p(class = "text-muted mb-0", "No condition added yet."))
      }
      tagList(lapply(seq_along(cond_state$ids), function(k) {
        cid       <- cond_state$ids[k]
        name_id   <- paste0("cond_name_", cid)
        col_id    <- paste0("cond_col_", cid)
        mod_id    <- paste0("cond_modalities_", cid)
        remove_id <- paste0("remove_cond_", cid)

        selected_col  <- isolate(input[[col_id]])
        selected_name <- isolate(input[[name_id]])
        if (is.null(selected_name) || !nzchar(selected_name)) selected_name <- paste0("Condition", k)
        modality_choices <- if (!is.null(selected_col) && nzchar(selected_col)) {
          clinic_modality_choices(selected_col)
        } else character(0)
        selected_mods <- isolate(input[[mod_id]])
        selected_mods <- selected_mods[selected_mods %in% modality_choices]

        div(
          class = "border rounded p-3 mb-2",
          tags$strong(paste("Condition", k)),
          textInput(ns(name_id), "Name :", value = selected_name),
          selectInput(ns(col_id), "Column :", choices = clinic_cols, selected = selected_col,
                      multiple = FALSE, selectize = FALSE, size = 5),
          selectInput(ns(mod_id), "Modalities :", choices = modality_choices,
                      selected = selected_mods, multiple = TRUE, selectize = FALSE, size = 5),
          div(class = "d-flex justify-content-end",
              actionButton(ns(remove_id), label = NULL, icon = icon("trash"),
                           class = "btn btn-sm btn-outline-danger"))
        )
      }))
    })

    register_condition <- function(cid) {
      local({
        current_cid <- cid
        current_col <- paste0("cond_col_", current_cid)
        current_mod <- paste0("cond_modalities_", current_cid)
        current_rm  <- paste0("remove_cond_", current_cid)
        observeEvent(input[[current_col]], {
          selected_col <- input[[current_col]]
          modality_choices <- if (!is.null(selected_col) && nzchar(selected_col)) {
            clinic_modality_choices(selected_col)
          } else character(0)
          selected_mods <- input[[current_mod]]
          selected_mods <- selected_mods[selected_mods %in% modality_choices]
          updateSelectInput(session, current_mod, choices = modality_choices,
                            selected = selected_mods)
        }, ignoreInit = TRUE)
        observeEvent(input[[current_rm]], {
          if (length(cond_state$ids) > 1) {
            cond_state$ids <- setdiff(cond_state$ids, current_cid)
          }
        }, ignoreInit = TRUE)
      })
    }

    add_condition <- function() {
      cid <- cond_state$next_id + 1
      cond_state$next_id <- cid
      cond_state$ids <- c(cond_state$ids, cid)
      register_condition(cid)
    }

    observeEvent(input$add_condition, { add_condition() })

    # Seed one condition by default (isolate: runs once, outside a reactive context).
    isolate(add_condition())

    conditions <- reactive({
      out <- list()
      for (k in seq_along(cond_state$ids)) {
        cid  <- cond_state$ids[k]
        name <- input[[paste0("cond_name_", cid)]]
        name <- trimws(if (is.null(name)) "" else as.character(name))
        if (!nzchar(name)) name <- paste0("Condition", k)
        col  <- input[[paste0("cond_col_", cid)]]
        col  <- trimws(if (is.null(col)) "" else as.character(col))
        mods <- input[[paste0("cond_modalities_", cid)]]
        mods <- trimws(as.character(mods))
        mods <- mods[!is.na(mods) & nzchar(mods)]
        if (!nzchar(col) || length(mods) == 0) next
        out[[length(out) + 1]] <- list(name = name, col = col, modalities = mods)
      }
      if (length(out) == 0) return(NULL)
      out
    })

    observeEvent(input$response_col, {
      req(input$response_col)
      mods <- clinic_modality_choices(input$response_col)
      updateSelectInput(session, "responders",     choices = mods, selected = character(0))
      updateSelectInput(session, "non_responders", choices = mods, selected = character(0))
    }, ignoreInit = TRUE)

    ## ---- Signatures: dynamic list + individual genes ----
    app_dir <- normalizePath(getwd())
    signatures_json <- fromJSON(file.path(app_dir, "REF_DATA", "signatures", "signatures.json"))

    observeEvent(input$bulk_file, {
      req(bulk_df())
      updateSelectizeInput(session, "genes", choices = rownames(bulk_df()),
                           selected = character(0), server = TRUE)
    })

    sig_state <- reactiveValues(ids = integer(), next_id = 0)

    output$signatures_ui <- renderUI({
      if (length(sig_state$ids) == 0) {
        return(tags$p(class = "text-muted mb-0", "No signature added yet."))
      }
      tagList(lapply(sig_state$ids, function(sid) {
        therapy_id   <- paste0("sig_therapy_", sid)
        signature_id <- paste0("sig_name_", sid)
        remove_id    <- paste0("remove_sig_", sid)
        sel_th <- isolate(input[[therapy_id]])
        sig_choices <- if (!is.null(sel_th) && nzchar(sel_th)) names(signatures_json[[sel_th]]) else character(0)
        div(
          class = "border rounded p-3 mb-2",
          selectizeInput(ns(therapy_id), "Therapy :", choices = names(signatures_json),
                         selected = sel_th, options = list(placeholder = "Therapy…")),
          selectizeInput(ns(signature_id), "Signature :", choices = sig_choices,
                         selected = isolate(input[[signature_id]]),
                         options = list(placeholder = "Signature…")),
          div(class = "d-flex justify-content-end",
              actionButton(ns(remove_id), label = NULL, icon = icon("trash"),
                           class = "btn btn-sm btn-outline-danger"))
        )
      }))
    })

    observeEvent(input$add_signature, {
      sid <- sig_state$next_id + 1
      sig_state$next_id <- sid
      sig_state$ids <- c(sig_state$ids, sid)
      local({
        current_sid  <- sid
        therapy_id   <- paste0("sig_therapy_", current_sid)
        signature_id <- paste0("sig_name_", current_sid)
        remove_id    <- paste0("remove_sig_", current_sid)
        observeEvent(input[[therapy_id]], {
          req(input[[therapy_id]])
          updateSelectizeInput(session, signature_id,
                               choices = names(signatures_json[[input[[therapy_id]]]]),
                               selected = isolate(input[[signature_id]]), server = TRUE)
        }, ignoreInit = TRUE)
        observeEvent(input[[remove_id]], {
          sig_state$ids <- setdiff(sig_state$ids, current_sid)
        }, ignoreInit = TRUE)
      })
    })

    selected_signatures <- reactive({
      sigs <- list()
      for (sid in sig_state$ids) {
        th  <- input[[paste0("sig_therapy_", sid)]]
        nm  <- input[[paste0("sig_name_", sid)]]
        if (is.null(th) || is.null(nm) || !nzchar(th) || !nzchar(nm)) next
        gs <- signatures_json[[th]][[nm]]$geneset
        if (is.null(gs) || length(gs) == 0) next
        sigs[[nm]] <- gs
      }
      for (g in input$genes) {
        if (!nzchar(g)) next
        sigs[[g]] <- stats::setNames(1, g)
      }
      if (length(sigs) == 0) return(NULL)
      names(sigs) <- make.unique(names(sigs))
      sigs
    })

    ## ---- Output directory ----
    shinyFiles::shinyDirChoose(input, id = "output_dir", session = session, roots = roots,
                               filetypes = c("", "txt", "tsv", "csv"))
    output_dir_path <- reactive({
      req(input$output_dir)
      shinyFiles::parseDirPath(roots, input$output_dir)
    })

    ## ---- Run ----
    comparison_res <- eventReactive(input$run_sigcomp, {
      req(clinic_df(), bulk_df(), selected_signatures(),
          conditions(),
          input$response_col, input$responders, input$non_responders,
          output_dir_path())

      withProgress(message = "Comparing signatures...", value = 0.05, {
        out_dir <- file.path(path.expand(output_dir_path()), "Signatures_comparison")
        res <- signatures_comparison_pipe(
          rnafilt_counts          = bulk_df(),
          clinic_annot            = clinic_df(),
          output_dir              = out_dir,
          signatures              = selected_signatures(),
          conditions              = conditions(),
          response_col            = input$response_col,
          responder_modalities    = input$responders,
          nonresponder_modalities = input$non_responders,
          responder_label         = input$responder_label %||% "R",
          nonresponder_label      = input$nonresponder_label %||% "NR",
          clinic_filters          = clinic_filters(),
          corr_method             = input$corr_method,
          corr_fdr                = isTRUE(input$corr_fdr),
          do_vst                  = isTRUE(input$do_vst),
          clin_continuous         = input$clin_continuous,
          progress_cb             = function(frac, detail) setProgress(value = frac, detail = detail)
        )
        setProgress(value = 1, detail = "done")
        res
      })
    })

    ## ---- Outputs ----
    output$logs <- renderPrint({
      note <- clinic_note()
      if (!is.null(note)) cat(note, "\n\n")
      res <- tryCatch(comparison_res(), error = function(e) NULL)
      if (is.null(res)) {
        if (is.null(note)) cat("Run the analysis to populate the logs.\n")
      } else {
        cat("Signatures compared:", ncol(res$scores), "\n")
        cat("Correlation subsets:", length(res$corr), "\n")
        cat("Output:", res$output_path, "\n")
      }
    })

    output$boxplot <- renderImage({
      p <- comparison_res()$boxplot_path
      req(p, file.exists(p))
      list(src = p, contentType = "image/png", alt = "Signatures boxplots")
    }, deleteFile = FALSE)

    output$clinical_boxplot <- renderImage({
      p <- comparison_res()$clinical_boxplot_path
      req(p, file.exists(p))
      list(src = p, contentType = "image/png", alt = "Clinical variable boxplots")
    }, deleteFile = FALSE)

    output$roc <- renderImage({
      p <- comparison_res()$roc_path
      req(p, file.exists(p))
      list(src = p, contentType = "image/png", alt = "Signature ROC curves")
    }, deleteFile = FALSE)

    output$corr_tabs_ui <- renderUI({
      res <- comparison_res()
      req(res)
      do.call(tabsetPanel, lapply(seq_along(res$corr), function(i) {
        tabPanel(
          title = res$corr[[i]]$label,
          card(div(class = "responsive-plot-frame",
                   imageOutput(ns(paste0("corr_img_", i)), width = "100%", height = "100%")))
        )
      }))
    })

    # One correlation image renderer per returned subset. The subset count is
    # 3 * n_conditions + 3; register a generous, lazy pool of renderers (each
    # only evaluates when its tab image is actually shown).
    for (i in seq_len(60L)) {
      local({
        idx <- i
        output[[paste0("corr_img_", idx)]] <- renderImage({
          res <- comparison_res()
          req(res, length(res$corr) >= idx)
          p <- res$corr[[idx]]$path
          req(p, file.exists(p))
          list(src = p, contentType = "image/png", alt = paste("Correlation", idx))
        }, deleteFile = FALSE)
      })
    }

    output$stats_table <- DT::renderDT({
      df <- comparison_res()$stats
      req(df)
      df
    })
  })
}
