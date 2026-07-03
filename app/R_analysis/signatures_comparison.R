# Compares several signatures across two user-defined sample conditions:
#   - boxplots of every signature, split by condition + response, with Wilcoxon
#   - 9 correlation matrices between signatures over each condition/response subset
#   - ROC curves per signature and condition
# Scores are computed on VST-normalized bulk, like the Signature Projection module.

signatures_comparison_pipe <- function(rnafilt_counts, clinic_annot, output_dir,
                                        signatures,
                                        conditions,
                                        response_col, responder_modalities, nonresponder_modalities,
                                        responder_label = "R", nonresponder_label = "NR",
                                        clinic_filters = NULL,
                                        corr_method = "spearman", corr_fdr = FALSE,
                                        sample_ID_col = "ID_Patient", do_vst = TRUE,
                                        progress_cb = NULL) {
  library(reticulate)

  report <- function(frac, detail) {
    if (is.function(progress_cb)) progress_cb(frac, detail)
  }
  report(0.05, "preprocessing & clinic filtering")

  use_python("/opt/conda/envs/BulkTools/bin/python", required = TRUE)
  source_python("py/py_plots.py")

  ## ---- Validation ----
  if (is.null(signatures) || length(signatures) == 0) {
    stop("Select at least one signature or gene to compare.")
  }
  for (col in c(cond1_col, cond2_col, response_col)) {
    if (is.null(col) || length(col) != 1 || !nzchar(col)) {
      stop("Condition columns and the response column must each be a single selected column.")
    }
  }

  ## ---- Optional global clinic filtering ----
  filtered <- filtering_on_clinic_and_genes(
    clinic_annot, rnafilt_counts,
    clinic_filters   = clinic_filters,
    filter_by_gene   = NULL, keep_low_or_high = NULL, quantile_thr = NULL,
    parsed_design    = "signatures_comparison",
    output_dir       = output_dir,
    folder_name      = "signatures_comparison"
  )
  rnafilt_counts <- filtered$rnafilt_counts
  clinic_annot   <- ensure_clinic_sample_id_col(filtered$clinic_annot)
  output_path    <- filtered$output_path

  for (col in c(cond1_col, cond2_col, response_col)) {
    if (!col %in% colnames(clinic_annot)) {
      stop(sprintf("Column '%s' not found in clinic_annot.", col))
    }
  }

  ## ---- (optional) VST normalization + sample alignment ----
  report(0.3, if (isTRUE(do_vst)) "VST normalization" else "using bulk as-is (no VST)")
  if (isTRUE(do_vst)) {
    vst <- normVST_bulk(round(rnafilt_counts))
  } else {
    vst <- as.matrix(rnafilt_counts)
  }

  ids <- intersect(clinic_annot$ID_Patient, colnames(vst))
  if (length(ids) == 0) {
    stop("No sample shared between clinic_annot and the bulk matrix after filtering.")
  }
  clinic_annot <- clinic_annot[match(ids, clinic_annot$ID_Patient), , drop = FALSE]
  vst <- vst[, ids, drop = FALSE]

  ## ---- Compute one score column per signature ----
  report(0.5, sprintf("computing scores for %d signature(s)", length(signatures)))
  score_cols <- lapply(names(signatures), function(nm) {
    sc <- compute_signature_score(vst, signatures[[nm]], scale_by_sum_abs = TRUE)
    sc[ids, "score"]
  })
  names(score_cols) <- names(signatures)
  scores_df <- as.data.frame(score_cols, check.names = FALSE)
  rownames(scores_df) <- ids

  ## ---- Define conditions and response membership ----
  trim_chr <- function(x) trimws(as.character(x))

  cond1_vals <- trim_chr(clinic_annot[[cond1_col]])
  cond2_vals <- trim_chr(clinic_annot[[cond2_col]])
  resp_vals  <- trim_chr(clinic_annot[[response_col]])

  cond1_ids <- ids[cond1_vals %in% trim_chr(cond1_modalities)]
  cond2_ids <- ids[cond2_vals %in% trim_chr(cond2_modalities)]

  if (length(cond1_ids) == 0) stop(sprintf("Condition 1 (%s) matched no sample.", cond1_name))
  if (length(cond2_ids) == 0) stop(sprintf("Condition 2 (%s) matched no sample.", cond2_name))

  resp_ids    <- ids[resp_vals %in% trim_chr(responder_modalities)]
  nonresp_ids <- ids[resp_vals %in% trim_chr(nonresponder_modalities)]

  map_response <- function(sample_ids) {
    out <- rep(NA_character_, length(sample_ids))
    out[sample_ids %in% resp_ids]    <- responder_label
    out[sample_ids %in% nonresp_ids] <- nonresponder_label
    out
  }

  cond1_scores <- scores_df[cond1_ids, , drop = FALSE]
  cond2_scores <- scores_df[cond2_ids, , drop = FALSE]
  meta1 <- data.frame(response = map_response(cond1_ids), row.names = cond1_ids,
                      stringsAsFactors = FALSE)
  meta2 <- data.frame(response = map_response(cond2_ids), row.names = cond2_ids,
                      stringsAsFactors = FALSE)

  if (sum(!is.na(meta1$response)) == 0 && sum(!is.na(meta2$response)) == 0) {
    stop("No sample carries a responder / non-responder label after mapping; check the response column and modalities.")
  }

  ## ---- Run the Python panels ----
  report(0.7, "building boxplots, correlation matrices and ROC curves")
  res <- run_signatures_comparison(
    cond1_scores      = cond1_scores,
    cond2_scores      = cond2_scores,
    meta1             = meta1,
    meta2             = meta2,
    response_col      = "response",
    cond1_name        = cond1_name,
    cond2_name        = cond2_name,
    responder_label   = responder_label,
    nonresponder_label = nonresponder_label,
    corr_method       = corr_method,
    corr_fdr          = corr_fdr,
    out_dir           = output_path
  )

  report(0.95, "finalizing")

  # Normalize the per-subset correlation entries into a tidy R list.
  corr_list <- lapply(res$corr, function(entry) {
    list(label = as.character(entry$label),
         path  = as.character(entry$path),
         n     = as.integer(entry$n))
  })

  list(
    boxplot_path = as.character(res$boxplot),
    roc_path     = as.character(res$roc),
    corr         = corr_list,
    stats        = as.data.frame(res$stats),
    scores       = scores_df,
    output_path  = output_path
  )
}
