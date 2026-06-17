# shiny_modules/mod_Tximport.R

mod_tximport_ui <- function(id) {
  ns <- NS(id)

  page_sidebar(
    sidebar = sidebar(
      accordion(
        id = ns("tximport_acc"),
        open = FALSE,
        multiple = TRUE,

        accordion_panel(
          "Files & tximport parameters *",
          shinyFiles::shinyDirButton(ns("bulk_folder"), "Bulk folder : *", "Select folder"),
          shinyFiles::shinyDirButton(ns("output_dir"),  "Output directory : *", "Select folder"),
          textInput(ns("file_strip"),
                    "File suffix(es) to strip from sample name (separate with '|') :",
                    value = "_quantif.tsv|_quantif.txt|_quantif.csv"),
          selectInput(ns("type"), "tximport type :",
                      choices = c("salmon", "sailfish", "alevin", "kallisto", "rsem", "stringtie", "none"),
                      selected = "salmon"),
          checkboxInput(ns("txIn"), "txIn", value = TRUE),
          selectInput(ns("countsFromAbundance"), "countsFromAbundance :",
                      choices = c("no", "scaledTPM", "lengthScaledTPM", "dtuScaledTPM"),
                      selected = "lengthScaledTPM"),
          textInput(ns("geneIdCol"),    "geneIdCol :",    value = "Gene_name"),
          textInput(ns("txIdCol"),      "txIdCol :",      value = "Name"),
          textInput(ns("abundanceCol"), "abundanceCol :", value = "TPM"),
          textInput(ns("countsCol"),    "countsCol :",    value = "NumReads"),
          textInput(ns("lengthCol"),    "lengthCol :",    value = "Length"),
          checkboxInput(ns("ignoreTxVersion"), "ignoreTxVersion", value = FALSE)
        ),

        accordion_panel(
          "Filtering parameters",
          numericInput(ns("min_gene_count"),
                       "Min gene count (a gene must reach this in enough samples) :",
                       value = 15, min = 0, step = 1),
          numericInput(ns("min_n_samples"),
                       "Min fraction of samples reaching the min gene count :",
                       value = 0.33, min = 0, max = 1, step = 0.01),
          numericInput(ns("min_seq_depth"),
                       "Min sequencing depth (total reads) per sample :",
                       value = 30000000, min = 0, step = 1000000),
          checkboxGroupInput(ns("remove_gene_classes"),
                             "Remove gene classes before normalization :",
                             choices = c("RPS" = "RPS", "RPL" = "RPL", "MT" = "MT"),
                             selected = character(0))
        )
      ),
      actionButton(ns("run_tximport"), "Run merge & filter")
    ),

    navset_card_tab(
      nav_panel("Logs", card(verbatimTextOutput(ns("logs")))),
      nav_panel("Summary", card(verbatimTextOutput(ns("summary")))),
      nav_panel("Read depth per sample", card(div(class = "responsive-plot-frame", plotOutput(ns("depth_plot"), width = "100%", height = "100%"))))
    )
  )
}


mod_tximport_server <- function(id, roots = c(home = "~")) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    app_dir <- normalizePath(getwd())
    tx2gene_fallback_gtf <- file.path(app_dir, "REF_DATA", "gencode.v36.annotation.gtf.gz")

    shinyFiles::shinyDirChoose(
      input,
      id = "bulk_folder",
      session = session,
      roots = roots,
      filetypes = c("")
    )
    shinyFiles::shinyDirChoose(
      input,
      id = "output_dir",
      session = session,
      roots = roots,
      filetypes = c("")
    )

    bulk_folder_path <- reactive({
      req(input$bulk_folder)
      shinyFiles::parseDirPath(roots, input$bulk_folder)
    })

    output_dir_path <- reactive({
      req(input$output_dir)
      shinyFiles::parseDirPath(roots, input$output_dir)
    })

    or_NULL <- function(x) {
      if (is.null(x) || length(x) == 0 || all(!nzchar(x))) NULL else x
    }

    tximport_res <- eventReactive(input$run_tximport, {
      req(bulk_folder_path(), output_dir_path())

      withProgress(message = "Merging & filtering bulks...", value = 0.05, {
        res <- tximport_merge_pipe(
          bulk_folder            = bulk_folder_path(),
          output_dir             = output_dir_path(),
          file_strip             = input$file_strip,
          tx_type                = input$type,
          tx_txIn                = input$txIn,
          tx_countsFromAbundance = input$countsFromAbundance,
          tx_geneIdCol           = input$geneIdCol,
          tx_txIdCol             = input$txIdCol,
          tx_abundanceCol        = input$abundanceCol,
          tx_countsCol           = input$countsCol,
          tx_lengthCol           = input$lengthCol,
          tx_ignoreTxVersion     = input$ignoreTxVersion,
          min_gene_count         = input$min_gene_count,
          min_n_samples          = input$min_n_samples,
          min_seq_depth          = input$min_seq_depth,
          remove_gene_classes    = or_NULL(input$remove_gene_classes),
          tx2gene_fallback_gtf   = tx2gene_fallback_gtf,
          progress_cb            = function(frac, detail) {
            setProgress(value = frac, detail = detail)
          }
        )
        setProgress(value = 1, detail = "done")
        res
      })
    })

    output$logs <- renderPrint({
      req(tximport_res())
      names(tximport_res())
    })

    output$summary <- renderText({
      req(tximport_res())
      tximport_res()$summary
    })

    output$depth_plot <- renderPlot({
      p <- tximport_res()$plots$depth
      req(p)
      print(p)
    })
  })
}
