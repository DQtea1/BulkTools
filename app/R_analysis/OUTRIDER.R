# Wraps the OUTRIDER pipeline (filterExpression -> controlForConfounders ->
# fit -> computePvalues / computeZscores -> results + plots).
# Vignette: https://bioconductor.org/packages/release/bioc/vignettes/OUTRIDER/inst/doc/OUTRIDER.pdf

OUTRIDER_pipe <- function(rnafilt_counts,
                          clinic_annot,
                          output_dir,
                          samples_to_exclude = NULL,
                          confounders        = NULL,
                          volcano_samples    = NULL,
                          plot_genes         = NULL,
                          iterations         = 3,
                          label_column     = "response") {
  library(OUTRIDER)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(org.Hs.eg.db)
  library(ggplot2)

  incProgress(0.15)

  rnafilt_counts <- round(as.matrix(rnafilt_counts))

  if (!"ID_Patient" %in% colnames(clinic_annot)) {
    clinic_annot$ID_Patient <- rownames(clinic_annot)
  }
  rownames(clinic_annot) <- clinic_annot$ID_Patient
  clinic_annot$sampleID  <- clinic_annot$ID_Patient

  common_samples <- intersect(rownames(clinic_annot), colnames(rnafilt_counts))
  if (length(common_samples) == 0) {
    stop("No samples shared between rnafilt_counts columns and clinic_annot.")
  }

  if (!is.null(confounders) && length(confounders) > 0) {
    missing_cols <- setdiff(confounders, colnames(clinic_annot))
    if (length(missing_cols) > 0) {
      stop(sprintf("Confounder column(s) not in clinic_annot: %s",
                   paste(missing_cols, collapse = ", ")))
    }
    keep <- complete.cases(clinic_annot[common_samples, confounders, drop = FALSE])
    common_samples <- common_samples[keep]
    if (length(common_samples) == 0) {
      stop("All samples were dropped because of NA values in the selected confounders.")
    }
  }

  clinic_annot   <- clinic_annot[common_samples, , drop = FALSE]
  rnafilt_counts <- rnafilt_counts[, common_samples, drop = FALSE]

  ods <- OutriderDataSet(countData = rnafilt_counts, colData = clinic_annot)


  ### PREPROCESSING ###
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  map <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys     = keys(txdb, keytype = "GENEID"),
    keytype  = "ENTREZID",
    columns  = c("SYMBOL")
  )

  ods <- filterExpression(ods, txdb, mapping = map, filterGenes = FALSE, savefpkm = TRUE)
  ods_filt <- ods[mcols(ods)$passedFilter, ]

  incProgress(0.20)

  ### Controlling for confounders ###
  ods_filt <- estimateSizeFactors(ods_filt)
  ods_filt <- estimateBestQ(ods_filt, useOHT = TRUE)
  ods_filt <- controlForConfounders(ods_filt, iterations = iterations)

  incProgress(0.30)

  ### Sample exclusion ###
  sampleExclusionMask(ods_filt) <- FALSE
  if (!is.null(samples_to_exclude) && length(samples_to_exclude) > 0) {
    samples_to_exclude <- intersect(samples_to_exclude, colnames(ods_filt))
    if (length(samples_to_exclude) > 0) {
      sampleExclusionMask(ods_filt[, samples_to_exclude]) <- TRUE
    }
  }


  ### Negative binomial fit + p-values + z-scores ###
  ods_filt <- fit(ods_filt)
  ods_filt <- computePvalues(ods_filt, alternative = "two.sided", method = "BY")
  ods_filt <- computeZscores(ods_filt)


  ### Results ###
  res <- results(ods_filt)
  res_df <- as.data.frame(res)

  tables_dir <- file.path(output_dir, "OUTRIDER", "tables")
  if (!dir.exists(tables_dir)) dir.create(tables_dir, recursive = TRUE)
  write.csv(res_df, file = file.path(tables_dir, "results_OUTRIDER_advanced.csv"), row.names = FALSE)

  incProgress(0.45)

  ### Aberrant per sample ###
  aberrant_plot <- plotAberrantPerSample(ods_filt, padjCutoff = 0.3)


  ### Volcano plots ###
  volcano_dir <- file.path(output_dir, "OUTRIDER", "plots", "volcanoes")
  if (!dir.exists(volcano_dir)) dir.create(volcano_dir, recursive = TRUE)

  if (is.null(volcano_samples) || length(volcano_samples) == 0) {
    volcano_samples <- unique(as.character(res_df$sampleID))
  }

  for (sample in volcano_samples) {
    p_volcano <- plotVolcano(ods_filt, sample, basePlot = TRUE) +
      labs(title = sample)
    ggsave(file.path(volcano_dir, paste0("volcano_", sample, ".png")), p_volcano)
  }

  incProgress(0.6)

  ### Gene-level plots ###
  rank_dir <- file.path(output_dir, "OUTRIDER", "plots", "expression_rank")
  exp_dir  <- file.path(output_dir, "OUTRIDER", "plots", "expected_vs_observed")
  if (!dir.exists(rank_dir)) dir.create(rank_dir, recursive = TRUE)
  if (!dir.exists(exp_dir))  dir.create(exp_dir,  recursive = TRUE)

  if (is.null(plot_genes) || length(plot_genes) == 0) {
    plot_genes <- unique(as.character(res_df$geneID))
  }

  ## Align clinic_annot to ods_filt sample order before labelling ##
  sample_order <- as.character(colnames(ods_filt))
  use_label <- !is.null(label_column) &&
               nzchar(label_column) &&
               label_column %in% colnames(clinic_annot)

  if (use_label) {
    clinic_aligned <- clinic_annot[match(sample_order, clinic_annot$ID_Patient), , drop = FALSE]
    label_vec <- clinic_aligned[[label_column]]
    colData(ods_filt)$label_col <- factor(
      label_vec,
      levels = unique(label_vec[!is.na(label_vec)])
    )
    groups_vec <- as.character(colData(ods_filt)$label_col)
  } else {
    groups_vec <- NULL
  }

  for (gene in plot_genes) {
    rank_args <- list(
      ods_filt,
      geneID   = gene,
      basePlot = TRUE,
      log      = TRUE,
      norm     = TRUE
    )
    exp_args <- list(
      ods_filt,
      geneID   = gene,
      basePlot = TRUE
    )
    if (!is.null(groups_vec)) {
      rank_args$groups <- groups_vec
      exp_args$groups  <- groups_vec
    }

    p_rank <- do.call(plotExpressionRank, rank_args)
    ggsave(file.path(rank_dir, paste0(gene, ".png")), p_rank)

    p_exp <- do.call(plotExpectedVsObservedCounts, exp_args)
    ggsave(file.path(exp_dir, paste0(gene, ".png")), p_exp)
  }

  incProgress(0.90)

  return(list(
    plots  = list(aberrant_per_sample = aberrant_plot),
    tables = list(results = res_df)
  ))
}
