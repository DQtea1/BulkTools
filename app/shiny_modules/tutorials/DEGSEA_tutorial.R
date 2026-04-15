degsea_tutorial_card <- function() {
  card(
    div(
      class = "p-3",
      tags$h3("How To Use DEGSEA"),
      tags$p(
        class = "mb-3",
        "The DEGSEA module lets you define two groups of bulk RNA-seq samples,",
        "compare them with DESeq2, and generate pathway-level summaries with",
        "GSEA, GESECA, and ssGSEA."
      ),

      tags$h4("What The Module Can Do"),
      tags$ul(
        tags$li("Load one bulk expression matrix and one clinical annotation table."),
        tags$li("Optionally pre-filter the dataset on one or more clinical variables."),
        tags$li("Optionally pre-filter the dataset on one gene-expression threshold."),
        tags$li("Build the control and test groups from multiple clinical filters and multiple gene-expression rules."),
        tags$li("Run differential expression and pathway analyses on the resulting sample selection.")
      ),

      tags$h4("Required Inputs And Where To Fill Them"),
      tags$ol(
        tags$li(
          tags$strong("Files *"),
          tags$ul(
            tags$li(tags$strong("Bulk path"), ": upload the bulk RNA-seq matrix."),
            tags$li(tags$strong("Clinic path"), ": upload the clinical annotation table."),
            tags$li(tags$strong("Output directory"), ": choose where all result folders and files will be written.")
          )
        ),
        tags$li(
          tags$strong("Filtering on clinic"),
          tags$ul(
            tags$li("Optional global pre-filter applied before DEGSEA."),
            tags$li("Use the ", tags$code("+"), " button to add one or more clinical filters."),
            tags$li("Within one row, multiple modalities are combined with OR."),
            tags$li("Across rows, filters are combined with AND.")
          )
        ),
        tags$li(
          tags$strong("Filtering on gene"),
          tags$ul(
            tags$li("Optional global pre-filter applied before group definition."),
            tags$li("Choose one gene, one quantile between 0 and 1, then keep either low or high expression samples.")
          )
        ),
        tags$li(
          tags$strong("Parameters *"),
          tags$ul(
            tags$li(tags$strong("DESeq covariates"), ": optional covariates added to the DESeq2 design."),
            tags$li(tags$strong("Control group definition *"), ": use clinical filters and optional expression filters to define the control samples."),
            tags$li(tags$strong("Test group definition *"), ": use clinical filters and optional expression filters to define the test samples."),
            tags$li("For control and test definitions, rows are combined by intersection."),
            tags$li(tags$strong("GSEA geneset *"), ": choose the pathway collection used for GSEA, GESECA, and ssGSEA.")
          )
        ),
        tags$li(
          tags$strong("Run DEGSEA"),
          tags$ul(
            tags$li("Click the button once all required inputs are filled.")
          )
        )
      ),

      tags$h4("Outputs In The App"),
      tags$ul(
        tags$li(tags$strong("Logs"), ": shows the returned objects available after the run."),
        tags$li(tags$strong("Volcano"), ": volcano plot of differential expression results."),
        tags$li(tags$strong("Differential Expression"), ": top differentially expressed genes shown as a table."),
        tags$li(tags$strong("GSEA table"), ": pathway enrichment scores and leading-edge genes."),
        tags$li(tags$strong("GESECA table"), ": pathway variance-based results."),
        tags$li(tags$strong("ssGSEA table"), ": per-sample enrichment score matrix.")
      ),

      tags$h4("Saved Results And File Formats"),
      tags$p(
        class = "mb-2",
        "Results are saved inside the selected output directory, in analysis-specific subfolders such as ",
        tags$code("DESeq/"),
        ", ",
        tags$code("GSEA/"),
        ", ",
        tags$code("GESECA/"),
        ", and ",
        tags$code("ssGSEA/"),
        "."
      ),
      tags$ul(
        tags$li(tags$code("DESeq_dds.rds"), ": compressed R object saved with ", tags$code("saveRDS(..., compress = \"xz\")"), "."),
        tags$li(tags$code("top_genes_down_500.csv"), " and ", tags$code("top_genes_up_500.csv"), ": CSV tables of top DE genes."),
        tags$li(tags$code("volcano_plot.jpg"), ": JPG image of the volcano plot."),
        tags$li(tags$code("ranked_genes.rnk"), ": tabular ranked-gene file used for pathway enrichment."),
        tags$li(tags$code("es_<geneset>.csv"), ": CSV file containing GSEA results."),
        tags$li(tags$code("<geneset>.csv"), ": CSV file containing GESECA results."),
        tags$li(tags$code("es_ssgsea_<geneset>.csv"), ": CSV file containing the ssGSEA enrichment matrix.")
      ),

      tags$h4("Practical Notes"),
      tags$ul(
        tags$li("A sample cannot belong to both control and test at the same time."),
        tags$li("If no sample matches a control or test definition, the run stops with an explicit error."),
        tags$li("Gene-expression rules in the control/test design are optional and are intersected with the clinical rules.")
      )
    )
  )
}
