ICA_tutorial_card <- function() {
  card(
    div(
      class = "p-3",
      tags$h3("How To Use ICA"),
      tags$p(
        class = "mb-3",
        "This module decomposes a bulk expression matrix with stabilized Independent",
        "Component Analysis. Each component gives a metagene (gene weights, matrix ",
        tags$code("S"), ") and its activity per sample (matrix ", tags$code("A"),
        "). The module has two steps: first estimate how many components to keep,",
        "then run the ICA and explore the results."
      ),

      tags$h4("Step 1 - Most Stable Transcriptomic Dimension (MSTD)"),
      tags$ul(
        tags$li("Load the ", tags$strong("Expression matrix"), " (genes x samples) and pick an ", tags$strong("Output directory"), "."),
        tags$li("Set the range of dimensions to scan (", tags$strong("Min / Max / Step"), ") and the number of ICA runs per dimension."),
        tags$li("Click ", tags$strong("Run MSTD"), ". The ", tags$strong("MSTD"), " tab shows the stability of the components as a function of the number of dimensions; the elbow / most stable region indicates a good number of components.")
      ),

      tags$h4("Step 2 - Run ICA"),
      tags$ul(
        tags$li("Set ", tags$strong("Number of components"), " to the value chosen from the MSTD plot, plus the number of runs and the FastICA algorithm (parallel / deflation)."),
        tags$li("Choose the two components for the ", tags$strong("Scatter with marginals"), " plot, and optionally a clinical column to color the points (requires a clinic file)."),
        tags$li("Click ", tags$strong("Run ICA"), ". The A (activities) and S (metagenes) matrices are saved under ", tags$code("output_dir/ICA/"), " as CSV + RDS.")
      ),

      tags$h4("Inputs"),
      tags$ul(
        tags$li(tags$strong("Expression matrix *"), ": counts / expression, genes in rows, samples in columns."),
        tags$li(tags$strong("Clinic (optional)"), ": sample annotation (row names = sample IDs, or an ", tags$code("ID_Patient"), " column). Enables scatter coloring and the clinical-association heatmaps."),
        tags$li(tags$strong("VST-normalize"), ": apply ", tags$code("normVST_bulk"), " before ICA (recommended for raw counts).")
      ),

      tags$h4("Output tabs"),
      tags$ul(
        tags$li(tags$strong("MSTD"), ": stability vs number of components (step 1)."),
        tags$li(tags$strong("Scatter with marginals"), ": two component activities with marginal histograms, colored by a clinical variable when provided."),
        tags$li(tags$strong("A matrix heatmap"), ": component activity per sample."),
        tags$li(tags$strong("Source distributions"), ": weight distribution and top genes of each metagene."),
        tags$li(tags$strong("Stability index"), ": reproducibility of each component."),
        tags$li(tags$strong("Component correlation"), ": correlation between component activities."),
        tags$li(tags$strong("Clinical association"), ": one heatmap per variable type - Spearman correlation for continuous variables, and -log10(p) (Mann-Whitney / Kruskal-Wallis) for categorical variables, with significance stars."),
        tags$li(tags$strong("Matrices"), ": the A and S tables (also saved to disk).")
      ),

      tags$h4("Practical Notes"),
      tags$ul(
        tags$li("Changing the scatter components or the color variable is applied on the next ", tags$strong("Run ICA"), "."),
        tags$li("A high number of runs improves stability estimates but is slower."),
        tags$li("Clinical association only appears when a clinic file with matching sample IDs is provided.")
      )
    )
  )
}
