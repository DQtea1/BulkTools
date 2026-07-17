signature_projection_tutorial_card <- function() {
  card(
    div(
      class = "p-3",
      tags$h3("How To Use Signature Projection"),
      tags$p(
        class = "mb-3",
        "This module evaluates a gene signature on a reference cohort and",
        "projects one new sample onto it. It learns, by nested cross-validation,",
        "where a responder / non-responder boundary sits on the signature score,",
        "then places the projected sample relative to that boundary with a",
        "calibrated probability of response and confidence interval."
      ),

      tags$h4("What The Module Can Do"),
      tags$ul(
        tags$li("Score a reference bulk cohort with a signature from ", tags$code("signatures.json"), " (or a combination of signatures / genes)."),
        tags$li("Project one external sample (from its Salmon quantification) onto the same score scale."),
        tags$li("Estimate the decision threshold and a gray zone by nested cross-validation, with bootstrap confidence intervals."),
        tags$li("Show a calibrated probability of response for the projected sample."),
        tags$li("Assess the performance of the signature with a ROC curve, a responder vs non-responder boxplot, and an optional Kaplan-Meier survival plot."),
        tags$li("Several signatures can be combined to refine the prediction accuracy.")

      ),

      tags$h4("Required Inputs And Where To Fill Them"),
      tags$ol(
        tags$li(
          tags$strong("Files *"),
          tags$ul(
            tags$li(tags$strong("Reference bulk"), ": bulk RNA-seq count matrix (genes x samples) used as the reference cohort."),
            tags$li(tags$strong("Reference clinic"), ": clinical annotation table for the reference samples. Sample IDs are read from the ", tags$code("ID_Patient"), " column when present, otherwise from row names."),
            tags$li(tags$strong("Sample to project"), ": the Salmon quantification file of the new sample to position; it is imported and merged with the reference bulk before normalization."),
            tags$li(tags$strong("Output directory"), ": where the figures are written.")
          )
        ),
        tags$li(
          tags$strong("Filtering on clinic"),
          tags$ul(
            tags$li("Optional pre-filter on the reference cohort. Use the ", tags$code("+"), " button to add clinical filters (multiple modalities in a row are OR, multiple rows are AND).")
          )
        ),
        tags$li(
          tags$strong("Filtering on gene"),
          tags$ul(
            tags$li("Optional: keep only the reference samples whose chosen gene expression sits below / above a given expression quantile (set quantile threshold to 0 disables it).")
          )
        ),
        tags$li(
          tags$strong("Select signature *"),
          tags$ul(
            tags$li(tags$strong("Signature family"), ": Select a signature family. (e.g. Ayers anti-PD1, CAF, Other immune...)"),
            tags$li(tags$strong("Signature"), ": Select a signature to use within the chosen family."),
            tags$li(tags$strong("Combine with other signatures / genes"), ": optionally add components combined with the base score using ", tags$code("+"), ", ", tags$code("-"), ", ", tags$code("*"), " or ", tags$code("/"), ". Each component is a signature or a single gene. Operators are applied left to right, in the order shown (no operation priority).")
          )
        ),
        tags$li(
          tags$strong("Parameters *"),
          tags$ul(
            tags$li(tags$strong("VST-normalize"), ": Whether to VST normalize the data. (only do this if the reference bulk contains raw counts)"),
            tags$li(tags$strong("Response Column"), ": a single clinical column holding the outcome."),
            tags$li(tags$strong("Responders / Non-Responders"), ": Select which modalities of the Response Column to include in each condition.")
          )
        ),
        tags$li(
          tags$strong("KM plot parameters"),
          tags$ul(
            tags$li(tags$strong("Plot KM plot ?"), ": set to ", tags$code("Ja"), " to also draw a Kaplan-Meier plot (default ", tags$code("Nein"), ", disabled)."),
            tags$li(tags$strong("Event / Survival time columns"), ": required only when the KM plot is enabled. Event realization column should contain 0/1 values (0 = censored, 1 = deceased), and Survival time column should contain numeric values."),
            tags$li(tags$strong("Quantile stratification"), ": median, tertile or quartile split of the score for the KM groups.")
          )
        ),
        tags$li(
          tags$strong("Run signature projection"),
          tags$ul(
            tags$li("Click once the files, the signature, the contrast and the responder / non-responder modalities are filled.")
          )
        )
      ),

      tags$h4("Outputs In The App"),
      tags$ul(
        tags$li(tags$strong("Logs"), ": shows the returned objects available after the run."),
        tags$li(tags$strong("Projection"), ": plots the reference score distribution, the decision threshold, the gray zone (uncertain classification), and the projected sample with its calibrated probability of response."),
        tags$li(tags$strong("Survival analysis"), ": plots the Kaplan-Meier plot (only when KM is enabled)."),
        tags$li(tags$strong("Signature assessment"), ": Evaluates the performance of the selected signature in discriminating responder and non-responder using boxplots (+Wilcoxon) and ROC curves with AUC.")
      ),

      tags$h4("Saved Results And File Formats"),
      tags$p(class = "mb-2", "Results are written under the output directory, in ", tags$code("<signature_family>/<signature>/"), " :"),
      tags$ul(
        tags$li(tags$code("projection/proj_<sample>.png"), ": the projection + confidence plot."),
        tags$li(tags$code("model_eval/ROC_AUC.png"), " and ", tags$code("model_eval/wilcox_boxplot.png"), ": the signature assessment plots."),
        tags$li(tags$code("KM_plot_<signature>_<quantile>.png"), ": the Kaplan-Meier plot, when enabled.")
      ),

      tags$h4("Practical Notes"),
      tags$ul(
        tags$li("Scores are computed on a VST normalization, as a weighted mean of the signature genes, so they are comparable with the Signatures Comparison module."),
        tags$li("The contrast must be a single column; the responder and non-responder selectors accept several modalities each."),
        tags$li("If a class has fewer than 2 samples after filtering, the run stops with a diagnostic explaining which step emptied it."),
        tags$li("The projected sample must share genes with the reference (same identifier type); a Salmon file with an empty gene column cannot be projected."),
        tags$li("Combining signatures works at the score level (the per-sample scores are combined), not at the gene level.")
      )
    )
  )
}