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
        tags$li("Set ", tags$strong("Number of components"), " to the value chosen from the MSTD plot, plus the number of stabilization runs."),
        tags$li("Choose the two metagenes for the ", tags$strong("Scatter with marginals"), " plot, and optionally a categorical clinical column to color the points (requires a clinic file)."),
        tags$li("Click ", tags$strong("Run ICA"), ". Two matrices are saved under ", tags$code("output_dir/ICA/"), " as CSV + RDS, following the notebook convention: ", tags$code("A_matrix_<K>metagenes"), " (gene weights, metagenes x genes) and ", tags$code("S_matrix_<K>metagenes"), " (sample activities, metagenes x samples).")
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
        tags$li(tags$strong("Scatter with marginals"), ": two metagene activities with stacked marginal histograms, colored by a categorical clinical variable when provided. The metagenes, colour column and bin count can be changed directly in the tab and update live (no ICA re-run)."),
        tags$li(tags$strong("Activity heatmap"), ": metagene activity per sample."),
        tags$li(tags$strong("Metagene gene weights"), ": gene-weight distribution and top genes of each metagene."),
        tags$li(tags$strong("Metagene correlation"), ": correlation between metagene activities."),
        tags$li(tags$strong("Clinical association"), ": one heatmap per variable type - Spearman correlation for continuous variables, and -log10(p) (Mann-Whitney / Kruskal-Wallis) for categorical variables, with significance stars."),
        tags$li(tags$strong("Matrices"), ": the A (gene weights) and S (sample activities) tables (also saved to disk).")
      ),

      tags$h4("Practical Notes"),
      tags$ul(
        tags$li("The scatter tab is interactive: changing the metagenes, colour or bins re-plots instantly from the fitted result, without re-running ICA."),
        tags$li("A high number of runs improves stability estimates but is slower."),
        tags$li("Clinical association only appears when a clinic file with matching sample IDs is provided.")
      )
    )
  )
}
