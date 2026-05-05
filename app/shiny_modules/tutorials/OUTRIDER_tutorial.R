outrider_tutorial_card <- function() {
  card(
    div(
      class = "p-3",
      tags$h3("How To Use OUTRIDER"),
      tags$p(
        class = "mb-3",
        tags$strong("Reference paper: "),
        tags$a(
          href = "https://pmc.ncbi.nlm.nih.gov/articles/PMC6288422/",
          target = "_blank",
          rel = "noopener noreferrer",
          "OUTRIDER: A Statistical Method for Detecting Aberrantly Expressed Genes in RNA Sequencing Data"
        )
      ),
      tags$p(
        class = "mb-3",
        "The OUTRIDER module flags genes whose expression is aberrant in",
        "individual samples. It runs the OUTRIDER pipeline (filter -> control",
        "for confounders -> negative binomial fit -> p-values + z-scores ->",
        "results + plots) and exposes the saved figures as browsable tabs."
      ),

      tags$h4("What The Module Can Do"),
      tags$ul(
        tags$li("Load one bulk count matrix and one clinical annotation table."),
        tags$li("Optionally mask one or more reference samples before fitting."),
        tags$li("Optionally control for known confounders during the autoencoder step."),
        tags$li("Pick a clinical column whose modalities will label samples in plots."),
        tags$li("Save volcano plots per sample and rank / expected-vs-observed plots per gene."),
        tags$li("Browse any saved figure from inside the app and open it as a sub-tab.")
      ),

      tags$h4("Required Inputs And Where To Fill Them"),
      tags$ol(
        tags$li(
          tags$strong("Files *"),
          tags$ul(
            tags$li(tags$strong("Bulk path"), ": upload the bulk RNA-seq count matrix (genes x samples). Counts must be integers; the module rounds them before fitting."),
            tags$li(tags$strong("Clinic path"), ": upload the clinical annotation table. Sample IDs are read from the ", tags$code("ID_Patient"), " column when present, otherwise from row names."),
            tags$li(tags$strong("Output directory"), ": choose where the result tables and figures will be written.")
          )
        ),
        tags$li(
          tags$strong("Parameters *"),
          tags$ul(
            tags$li(tags$strong("Samples to exclude"), ": List of sample IDs to mask before fitting. Useful when one sample is a known outlier or when a sample has replicates and you do not want those to influence the autoencoder."),
            tags$li(tags$strong("Confounders"), ": clinical columns whose effect should be regressed out by the autoencoder. Samples with ", tags$code("NA"), " in any selected confounder are dropped from the run."),
            tags$li(tags$strong("Labelizing column"), ": Clinical column whose modalities are used to color samples in the plots."),
            tags$li(tags$strong("Samples to volcano"), ": sample IDs for which a per-sample volcano plot is saved. Leave empty to generate one for every sample present in the results table."),
            tags$li(tags$strong("Genes to plot"), ": gene IDs for which an expression-rank plot and an expected-vs-observed plot are saved. Leave empty to generate them for every gene present in the results table (this can produce many files)."),
            tags$li(tags$strong("Iterations for confounder control"), ": number of iterations passed to ", tags$code("controlForConfounders"), ". Higher is more stable but slower.")
          )
        ),
        tags$li(
          tags$strong("Run OUTRIDER"),
          tags$ul(
            tags$li("Click the button once the bulk file, clinic file and output directory are filled.")
          )
        )
      ),

      tags$h4("Outputs In The App"),
      tags$ul(
        tags$li(tags$strong("Logs"), ": shows the returned objects available after the run."),
        tags$li(tags$strong("Aberrant per sample"), ": bar plot of the number of aberrantly expressed genes per sample."),
        tags$li(tags$strong("Results"), ": full OUTRIDER results table (sample, gene, p-value, padj, z-score, ...)."),
        tags$li(tags$strong("Volcano"), ": click the picker to browse ", tags$code("plots/volcanoes/"), " and open any saved volcano plot as a sub-tab."),
        tags$li(tags$strong("Expression rank"), ": click the picker to browse ", tags$code("plots/expression_rank/"), " and open any saved rank plot as a sub-tab."),
        tags$li(tags$strong("Expected vs observed"), ": click the picker to browse ", tags$code("plots/expected_vs_observed/"), " and open any saved scatter plot as a sub-tab.")
      ),
      tags$p(
        class = "text-muted small",
        "Each picked plot opens as its own sub-tab; selecting the same file again refocuses",
        "the existing sub-tab instead of opening a duplicate."
      ),

      tags$h4("Saved Results And File Formats"),
      tags$p(
        class = "mb-2",
        "Results are saved inside the selected output directory:"
      ),
      tags$ul(
        tags$li(tags$code("tables/results_OUTRIDER_advanced.csv"), ": full results table from ", tags$code("OUTRIDER::results()"), "."),
        tags$li(tags$code("plots/volcanoes/volcano_<sample>.png"), ": per-sample volcano plot."),
        tags$li(tags$code("plots/expression_rank/<gene>.png"), ": expression-rank plot for one gene."),
        tags$li(tags$code("plots/expected_vs_observed/<gene>.png"), ": expected vs observed counts for one gene.")
      ),

      tags$h4("Practical Notes"),
      tags$ul(
        tags$li("Sample IDs must overlap between the bulk count matrix columns and the clinic table; samples that are not in both are dropped."),
        tags$li("Selecting confounders with missing values implicitly removes the affected samples from the run."),
        tags$li("Leaving ", tags$strong("Samples to volcano"), " or ", tags$strong("Genes to plot"), " empty can produce a large number of figures; use the selectors to restrict to the sample/genes you actually want to inspect."),
        tags$li("The aberrant-per-sample plot uses ", tags$code("padjCutoff = 0.3"), "; tweak the OUTRIDER pipeline file if you need a different threshold."),
        tags$li("Plots are generated only at run time; the picker tabs only show files once OUTRIDER has finished running.")
      )
    )
  )
}
