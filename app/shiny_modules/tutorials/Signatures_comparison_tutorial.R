signatures_comparison_tutorial_card <- function() {
  card(
    div(
      class = "p-3",
      tags$h3("How To Use Signatures Comparison"),
      tags$p(
        class = "mb-3",
        "This module evaluates and compares several gene signatures (and/or",
        "individual genes) side by side. You define two sample groups",
        "(\"conditions\") and a responder / non-responder split, then the module",
        "scores every signature per sample and produces boxplots, correlation",
        "matrices and ROC curves to compare their behaviour."
      ),

      tags$h4("What The Module Can Do"),
      tags$ul(
        tags$li("Score any number of signatures and single genes on the same bulk cohort."),
        tags$li("Compare score distributions between responders and non-responders, per condition, with a Wilcoxon (Mann-Whitney) test."),
        tags$li("Inspect how signatures correlate with each other across 9 sample subsets."),
        tags$li("Read each signature's discriminative power (ROC / AUC) within each condition.")
      ),

      tags$h4("Required Inputs And Where To Fill Them"),
      tags$ol(
        tags$li(
          tags$strong("Files *"),
          tags$ul(
            tags$li(tags$strong("Bulk path"), ": bulk RNA-seq count matrix (genes x samples). Scores are computed on a VST normalization of these counts."),
            tags$li(tags$strong("Clinic path"), ": clinical annotation table. Sample IDs are read from the ", tags$code("ID_Patient"), " column when present, otherwise from row names."),
            tags$li(tags$strong("Output directory"), ": where the figures and tables are written.")
          )
        ),
        tags$li(
          tags$strong("Filter on clinic & define conditions"),
          tags$ul(
            tags$li(tags$strong("Global clinical filters"), ": optional pre-filter applied to every sample before anything else (use the ", tags$code("+"), " button; multiple modalities in a row are OR, multiple rows are AND)."),
            tags$li(tags$strong("Conditions"), ": one condition by default; click the ", tags$code("+"), " button to add as many as you want (trash icon to remove; at least one is kept). Each condition is a sample group: give it a name, pick a clinical column, then pick the modalities that define the group. Conditions are independent (e.g. cohorts, treatments, pre vs post)."),
            tags$li(tags$strong("Response column"), ": the clinical column holding the outcome, then select which of its modalities count as ", tags$strong("responders"), " and which as ", tags$strong("non-responders"), ".")
          )
        ),
        tags$li(
          tags$strong("Select signatures *"),
          tags$ul(
            tags$li("Use the ", tags$code("+"), " button to add as many signatures as you want (each is chosen by therapy then signature name, from ", tags$code("signatures.json"), ")."),
            tags$li("Add individual genes in the gene selector; each gene is treated as a one-gene signature.")
          )
        ),
        tags$li(
          tags$strong("Parameters"),
          tags$ul(
            tags$li(tags$strong("Correlation test"), ": spearman (default), pearson or kendall, used for the correlation matrices."),
            tags$li(tags$strong("BH-corrected q-values"), ": when ticked, the significance stars use Benjamini-Hochberg q-values instead of raw p-values (recommended when comparing many signatures)."),
            tags$li(tags$strong("Responder / Non-responder labels"), ": the labels shown on the figures (default R / NR).")
          )
        ),
        tags$li(
          tags$strong("Run comparison"),
          tags$ul(
            tags$li("Click once the files, at least one condition, the response split and at least one signature are filled.")
          )
        )
      ),

      tags$h4("Outputs In The App"),
      tags$ul(
        tags$li(tags$strong("Boxplots"), ": every signature side by side, split by condition and response status, with the Wilcoxon p-value comparing responders vs non-responders."),
        tags$li(tags$strong("Correlations"), ": one sub-tab per sample subset, each a correlation matrix between signatures — for every condition its responders, non-responders and whole group, plus all responders, all non-responders and all samples (", tags$code("3 x n_conditions + 3"), " subsets). Stars: *** p<0.001, ** p<0.01, * p<0.05."),
        tags$li(tags$strong("ROC curves"), ": a grid (one row per condition, one column per signature) showing each signature's ROC and AUC for separating responders from non-responders."),
        tags$li(tags$strong("Wilcoxon stats"), ": the full table of Mann-Whitney statistics and p-values per signature and condition.")
      ),

      tags$h4("Saved Results And File Formats"),
      tags$p(class = "mb-2", "Results are written inside the output directory, under ", tags$code("Signatures_comparison/"), " :"),
      tags$ul(
        tags$li(tags$code("signatures_boxplots.png"), ": the boxplot panel."),
        tags$li(tags$code("correlations/corr_<i>_<subset>.png"), ": the 9 correlation matrices."),
        tags$li(tags$code("signature_rocs.png"), ": the ROC grid.")
      ),

      tags$h4("Practical Notes"),
      tags$ul(
        tags$li("Scores use the same VST + weighted-mean convention as the Signature Projection module, so values are comparable between the two modules."),
        tags$li("Conditions are independent selections: if they overlap, shared samples are counted in both, and the 'all samples / all responders' correlation subsets will include them twice. Keep the conditions disjoint to avoid this."),
        tags$li("ROC orientation is raw (higher score = predicted responder); a signature that is anti-correlated with response will show an AUC below 0.5 rather than being auto-flipped."),
        tags$li("A correlation matrix needs at least 3 samples and 2 signatures in the subset; sparser subsets are shown as 'Not enough data'.")
      )
    )
  )
}