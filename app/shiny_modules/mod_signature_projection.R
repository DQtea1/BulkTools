# shiny_modules/mod_signature_projection.R
library(jsonlite)

mod_signature_proj_ui <- function(id) {
  ns <- NS(id)
  tagList(
    page_sidebar(
      sidebar = sidebar(
        accordion(
          id = ns("sign_proj_acc"),
          open = FALSE,
          multiple = TRUE,

          accordion_panel(
            "Files *",
            fileInput(ns("bulk_file"),   "Reference bulk : *",   accept = c(".csv", ".tsv")),
            fileInput(ns("clinic_file"), "Reference clinic : *", accept = c(".csv", ".tsv")),
            fileInput(ns("proj_sample"), "Sample to project : *", accept = c(".csv", ".tsv")),
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
            )
            ,
            uiOutput(ns("clinic_filters_ui"))
          ),
          accordion_panel(
            "Filtering on gene",
            selectizeInput(
              ns("gene_filt_proj"), "Filter on gene :",
              choices = stats::setNames("", ""),
              selected = "",
              multiple = FALSE,
              options = list(placeholder = "Type a gene symbol…")
            ),
            numericInput(
              ns("quantile_proj"), "Quantile threshold :",
              value = 0, min = 0, max = 1, step = 0.05
            ),
            selectInput(
              ns("low_or_high_proj"), "Keep low or high expression :",
              choices = c("low", "high"),
              multiple = FALSE, selectize = FALSE, size = 2
            )
          ),
          accordion_panel(
            "Select signature *",
            selectizeInput(
              ns("therapy"), "Therapy :",
              choices = character(0),
              multiple = FALSE,
              options = list(placeholder = "Type a treatment (no spaces)…")
            ),
            selectizeInput(
              ns("signatures"), "Signature :",
              choices = character(0),
              multiple = FALSE,
              options = list(placeholder = "Type a signature name (no spaces)…")
            ),
            div(
              class = "d-flex justify-content-between align-items-center mb-2 mt-2",
              tags$span("Combine with other signatures / genes"),
              actionButton(
                ns("add_component"),
                label = NULL,
                icon = icon("plus"),
                class = "btn btn-sm btn-outline-primary"
              )
            ),
            uiOutput(ns("extra_components_ui"))
          ),
          accordion_panel(
            "Parameters *",
            selectInput(
              ns("contrast"), "Contrast :",
              choices = character(0),
              multiple = FALSE, selectize = FALSE, size = 7
            ),
            selectInput(
              ns("responders"), "Responders :",
              choices = character(0),
              multiple = TRUE, selectize = FALSE, size = 7
            ),
            selectInput(
              ns("non_responders"), "Non-Responders :",
              choices = character(0),
              multiple = TRUE, selectize = FALSE, size = 7
            )
          ),
          accordion_panel(
            "KM plot parameters ",
            selectInput(
              ns("plot_km_or_not"), "Plot KM plot ? :",
              choices = c("Ja", "Nein"),
              selected = "Nein",
              multiple = FALSE, selectize = FALSE, size = 2
            ),
            selectInput(
              ns("event_realization_col"), "Event realization column :",
              choices = character(0),
              multiple = FALSE, selectize = FALSE, size = 7
            ),
            selectInput(
              ns("survival_time_col"), "Survival time column :",
              choices = character(0),
              multiple = FALSE, selectize = FALSE, size = 7
            ),
            selectInput(
              ns("group_quantile"), "Quantile stratification :",
              choices = c("median", "tertile", "quartile"),
              multiple = FALSE, selectize = FALSE, size = 3
            ),
          )
        ),
        actionButton(ns("run_signature_proj"), "Run signature projection")
      ),
      navset_card_tab(
        nav_panel("Logs", card(verbatimTextOutput(ns("logs")))),
        nav_panel("Projection",
                  card(div(class = "responsive-plot-frame",
                           imageOutput(ns("projection"), width = "100%", height = "100%")))),
        nav_panel("Survival analysis",
                  card(div(class = "responsive-plot-frame",
                           imageOutput(ns("KM_plot"), width = "100%", height = "100%")))),
        nav_panel("Signature assessment",
                  layout_columns(
                    col_width = c(6, 6, 6),
                    card(div(class = "responsive-plot-frame",
                             imageOutput(ns("signature_assessment_boxplot"), width = "100%", height = "100%"))),
                    card(div(class = "responsive-plot-frame",
                             imageOutput(ns("signature_assessment_ROC"), width = "100%", height = "100%")))
                    # card(div(class = "responsive-plot-frame",
                    #          imageOutput(ns("signature_assessment_conf_mat"), width = "100%", height = "100%")))
                  )
        )
      )    
    )
  )  
}


mod_signature_proj_server <- function(id, roots = c(home = "~")) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    library(reticulate)
    Sys.unsetenv("RETICULATE_PYTHON")
    use_python("/opt/conda/envs/BulkTools/bin/python", required = TRUE)
    source_python("py/py_plots.py")

    clinic_df <- reactive({
      read_delim_auto(input$clinic_file)
    })

    bulk_df <- reactive({
      read_delim_auto(input$bulk_file)
    })

    gene_filter_choices <- function(genes = character(0)) {
      c(stats::setNames("", ""), stats::setNames(genes, genes))
    }

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

    observeEvent(input$clinic_file, {
      req(clinic_df())
      updateSelectInput(session, "contrast", choices = names(clinic_df()), selected = "")
      updateSelectInput(session, "survival_time_col", choices = names(clinic_df()),  selected = "") 
      updateSelectInput(session, "event_realization_col", choices = names(clinic_df()),  selected = "") 
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

    # update gene list when bulk changes
    observeEvent(input$bulk_file, {
      req(bulk_df())
      genes <- rownames(bulk_df())
      updateSelectizeInput(session, "gene_filt_proj", choices = gene_filter_choices(genes), selected = "", server = TRUE)
    })

    observeEvent(input$contrast, {
      req(clinic_df(), input$contrast)
      v <- clinic_df()[[input$contrast]]
      mods <- if (is.factor(v)) levels(v) else sort(unique(as.character(v)))
      mods <- mods[!is.na(mods) & mods != ""]
      updateSelectInput(session, "responders", choices = mods, selected = "")
      updateSelectInput(session, "non_responders",    choices = mods, selected = "")
    })

    app_dir <- normalizePath(getwd())
    signatures_path = file.path(app_dir, "REF_DATA", "signatures", "signatures.json")
    signatures_json = fromJSON(signatures_path)
    updateSelectizeInput(session, "therapy", choices = names(signatures_json), server = TRUE)

    observeEvent(input$therapy, {
      req(input$therapy, signatures_json)
      signatures_list <- names(signatures_json[[input$therapy]])
      updateSelectizeInput(session, "signatures", choices = signatures_list, server = TRUE)
    })

    sel_signature <- reactive({
      req(input$therapy, input$signatures, signatures_json)
      signatures_json[[input$therapy]][[input$signatures]]
    })

    ## ---- Extra signature / gene components (combine with +, -, *, /) ----
    `%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

    extra_comp_state <- reactiveValues(ids = integer(), next_id = 0)

    comp_gene_choices <- function() {
      if (isTruthy(input$bulk_file)) rownames(bulk_df()) else character(0)
    }

    output$extra_components_ui <- renderUI({
      if (length(extra_comp_state$ids) == 0) {
        return(tags$p(class = "text-muted mb-0", "No extra component."))
      }
      tagList(lapply(extra_comp_state$ids, function(cid) {
        op_id        <- paste0("comp_op_", cid)
        type_id      <- paste0("comp_type_", cid)
        therapy_id   <- paste0("comp_therapy_", cid)
        signature_id <- paste0("comp_signature_", cid)
        gene_id      <- paste0("comp_gene_", cid)
        remove_id    <- paste0("remove_component_", cid)

        div(
          class = "border rounded p-3 mb-2",
          selectInput(ns(op_id), "Operator :",
                      choices = c("add (+)" = "+", "subtract (-)" = "-",
                                  "multiply (*)" = "*", "divide (/)" = "/"),
                      selected = isolate(input[[op_id]]) %||% "-",
                      multiple = FALSE, selectize = FALSE, size = 4),
          selectInput(ns(type_id), "Component type :",
                      choices = c("Signature", "Gene"),
                      selected = isolate(input[[type_id]]) %||% "Signature",
                      multiple = FALSE, selectize = FALSE, size = 2),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'Signature'", ns(type_id)),
            selectizeInput(ns(therapy_id), "Therapy :",
                           choices = names(signatures_json),
                           selected = isolate(input[[therapy_id]]),
                           options = list(placeholder = "Therapy…")),
            selectizeInput(ns(signature_id), "Signature :",
                           choices = character(0),
                           selected = isolate(input[[signature_id]]),
                           options = list(placeholder = "Signature…"))
          ),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'Gene'", ns(type_id)),
            selectizeInput(ns(gene_id), "Gene :",
                           choices = character(0),
                           selected = isolate(input[[gene_id]]),
                           options = list(placeholder = "Type a gene symbol…"))
          ),
          div(class = "d-flex justify-content-end",
              actionButton(ns(remove_id), label = NULL, icon = icon("trash"),
                           class = "btn btn-sm btn-outline-danger"))
        )
      }))
    })

    observeEvent(input$add_component, {
      cid <- extra_comp_state$next_id + 1
      extra_comp_state$next_id <- cid
      extra_comp_state$ids <- c(extra_comp_state$ids, cid)

      local({
        current_cid  <- cid
        therapy_id   <- paste0("comp_therapy_", current_cid)
        signature_id <- paste0("comp_signature_", current_cid)
        remove_id    <- paste0("remove_component_", current_cid)

        observeEvent(input[[therapy_id]], {
          req(input[[therapy_id]])
          sigs <- names(signatures_json[[input[[therapy_id]]]])
          updateSelectizeInput(session, signature_id, choices = sigs,
                               selected = isolate(input[[signature_id]]), server = TRUE)
        }, ignoreInit = TRUE)

        observeEvent(input[[remove_id]], {
          extra_comp_state$ids <- setdiff(extra_comp_state$ids, current_cid)
        }, ignoreInit = TRUE)
      })
    })

    # Fill the per-component gene selectize inputs server-side (the gene list is
    # large, so client-side rendering is slow). Runs after each UI flush.
    observe({
      ids <- extra_comp_state$ids
      genes <- comp_gene_choices()
      session$onFlushed(function() {
        for (cid in ids) {
          gene_id <- paste0("comp_gene_", cid)
          updateSelectizeInput(
            session, gene_id,
            choices  = genes,
            selected = isolate(input[[gene_id]]) %||% "",
            server   = TRUE
          )
        }
      }, once = TRUE)
    })

    extra_components <- reactive({
      if (length(extra_comp_state$ids) == 0) return(NULL)
      comps <- list()
      for (cid in extra_comp_state$ids) {
        op   <- input[[paste0("comp_op_", cid)]]
        type <- input[[paste0("comp_type_", cid)]]
        if (is.null(op) || is.null(type)) next

        if (identical(type, "Signature")) {
          th  <- input[[paste0("comp_therapy_", cid)]]
          sig <- input[[paste0("comp_signature_", cid)]]
          if (is.null(th) || is.null(sig) || !nzchar(th) || !nzchar(sig)) next
          gs <- signatures_json[[th]][[sig]]$geneset
        } else {
          g <- input[[paste0("comp_gene_", cid)]]
          if (is.null(g) || !nzchar(g)) next
          gs <- stats::setNames(1, g)
        }
        if (is.null(gs) || length(gs) == 0) next
        comps[[length(comps) + 1]] <- list(op = op, geneset = gs)
      }
      if (length(comps) == 0) NULL else comps
    })

    projection_res <- eventReactive(input$run_signature_proj, {
      req(
        clinic_df(), bulk_df(), sel_signature(),
        input$proj_sample$datapath, output_dir_path()
      )

      withProgress(message = "Projecting signature score...", value = 0.05, {
        output_dir = file.path(path.expand(output_dir_path()), input$therapy, input$signatures)
        res <- signature_proj_pipe(
          rnafilt_counts         = bulk_df(),
          clinic_annot           = clinic_df(),
          output_dir             = output_dir,
          therapy_used           = input$therapy,
          signature_name         = input$signatures,
          signature_to_use       = sel_signature()$geneset,
          sample_to_project_path = input$proj_sample$datapath[1],
          contrast               = input$contrast,
          resp_var               = input$responders,
          non_resp_var           = input$non_responders,
          survival_time_col      = input$survival_time_col,
          event_realization_col  = input$event_realization_col,
          group_quantile         = input$group_quantile,
          do_km_plot             = identical(input$plot_km_or_not, "Ja"),
          clinic_filters         = clinic_filters(),
          filter_by_gene         = input$gene_filt_proj,
          keep_low_or_high       = input$low_or_high_proj,
          quantile_thr           = input$quantile_proj,
          extra_components       = extra_components(),
          progress_cb            = function(frac, detail) {
            setProgress(value = frac, detail = detail)
          }
        )
        setProgress(value = 1, detail = "done")
        res
      })
    })

    output$logs <- renderPrint({
      req(projection_res())
      names(projection_res())
    })

    output$projection <- renderImage({
      width  <- session$clientData$output_myplot_width
      height <- session$clientData$output_myplot_height
      p <- projection_res()$proj_plot$save_path
      req(p, file.exists(p))
      list(
        src = p,
        contentType = "image/png",
        alt = "Signature projection plot"
      )
    }, deleteFile = FALSE)

    output$KM_plot <- renderImage({
      width  <- session$clientData$output_myplot_width
      height <- session$clientData$output_myplot_height
      p <- projection_res()$KM_plot$save_path
      req(p, file.exists(p))
      list(
        src = p,
        width = width,
        height = height,
        contentType = "image/png",
        alt = "Signature survival KM"
      )
    }, deleteFile = FALSE)

    output$signature_assessment_ROC <- renderImage({
      width  <- session$clientData$output_myplot_width
      height <- session$clientData$output_myplot_height
      roc_p <- projection_res()$roc_plot$roc_save_path
      req(roc_p, file.exists(roc_p))
      list(
        src = roc_p,
        width = width,
        height = height,
        contentType = "image/png",
        alt = "Signature evaluation ROC"
      )

    }, deleteFile = FALSE)
    output$signature_assessment_boxplot <- renderImage({
      width  <- session$clientData$output_myplot_width
      height <- session$clientData$output_myplot_height
      box_p <- projection_res()$box_plot$save_path
      req(box_p, file.exists(box_p))
      list(
        src = box_p,
        width = width,
        height = height,
        contentType = "image/png",
        alt = "Signature evaluation boxplot"
      )
    }, deleteFile = FALSE)

  })
}
