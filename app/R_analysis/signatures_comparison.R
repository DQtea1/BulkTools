# Compares several signatures across any number of user-defined sample conditions:
#   - boxplots of every signature, split by condition + response, with Wilcoxon
#   - correlation matrices between signatures over each condition/response subset
#   - ROC curves per signature and condition
# Scores are computed on (optionally VST-normalized) bulk, like Signature Projection.
# `conditions` is a list of list(name=, col=, modalities=) entries.

signatures_comparison_pipe <- function(rnafilt_counts, clinic_annot, output_dir,
                                        signatures,
                                        conditions,
                                        response_col, responder_modalities, nonresponder_modalities,
                                        responder_label = "R", nonresponder_label = "NR",
                                        clinic_filters = NULL,
                                        corr_method = "spearman", corr_fdr = FALSE,
                                        sample_ID_col = "ID_Patient", do_vst = TRUE,
                                        clin_continuous = NULL,
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
  if (is.null(conditions) || length(conditions) == 0) {
    stop("Define at least one condition.")
  }
  if (is.null(response_col) || length(response_col) != 1 || !nzchar(response_col)) {
    stop("The response column must be a single selected column.")
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

  if (!response_col %in% colnames(clinic_annot)) {
    stop(sprintf("Column '%s' not found in clinic_annot.", response_col))
  }

  ## ---- (optional) VST normalization + sample alignment ----
  report(0.3, if (isTRUE(do_vst)) "VST normalization" else "using bulk as-is (no VST)")
  if (isTRUE(do_vst)) {
    vst <- normVST_bulk(round(rnafilt_counts))
  } else {
    vst <- as.matrix(rnafilt_counts)
  }

  # Trim IDs on both sides (whitespace is a common cause of a false mismatch) and
  # write them back so the subsequent match()/subset use the same trimmed values.
  clinic_annot$ID_Patient <- trimws(as.character(clinic_annot$ID_Patient))
  colnames(vst)           <- trimws(as.character(colnames(vst)))
  clinic_ids_all <- clinic_annot$ID_Patient
  bulk_ids_all   <- colnames(vst)
  ids <- intersect(clinic_ids_all, bulk_ids_all)
  if (length(ids) == 0) {
    stop(sprintf(
      paste0(
        "No sample shared between clinic_annot and the bulk matrix after filtering.\n",
        "--- DIAGNOSTIC ---\n",
        "clinic rows (after NA-ID + clinic filters) : %d\n",
        "bulk samples (columns)                     : %d\n",
        "clinic ID_Patient (first 5)                : %s\n",
        "bulk column names (first 5)                : %s\n",
        "------------------\n",
        "Common causes: (a) the clinic sample IDs are in a different column than the ",
        "bulk column names, (b) an ID-type mismatch (whitespace, case, separators, ",
        "prefix/suffix), (c) the global clinic filter or NA-ID drop removed every matching sample."
      ),
      length(clinic_ids_all),
      length(bulk_ids_all),
      paste(head(clinic_ids_all, 5), collapse = ", "),
      paste(head(bulk_ids_all, 5), collapse = ", ")
    ))
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

  ## ---- Response membership ----
  trim_chr <- function(x) trimws(as.character(x))

  resp_vals   <- trim_chr(clinic_annot[[response_col]])
  resp_ids    <- ids[resp_vals %in% trim_chr(responder_modalities)]
  nonresp_ids <- ids[resp_vals %in% trim_chr(nonresponder_modalities)]

  map_response <- function(sample_ids) {
    out <- rep(NA_character_, length(sample_ids))
    out[sample_ids %in% resp_ids]    <- responder_label
    out[sample_ids %in% nonresp_ids] <- nonresponder_label
    out
  }

  ## ---- Optional: numeric matrix of the selected continuous clinical variables ----
  cont_vars <- character(0)
  clin_num  <- NULL
  if (!is.null(clin_continuous) && length(clin_continuous) > 0) {
    cont_vars <- intersect(as.character(clin_continuous), colnames(clinic_annot))
    if (length(cont_vars) > 0) {
      clin_num <- clinic_annot[, cont_vars, drop = FALSE]
      clin_num[] <- lapply(clin_num, function(x) suppressWarnings(as.numeric(as.character(x))))
      rownames(clin_num) <- clinic_annot$ID_Patient
    }
  }

  ## ---- Build one score/meta subset per condition (any number) ----
  cond_names       <- character(0)
  cond_scores_list <- list()
  cond_meta_list   <- list()
  cond_ids_list    <- list()
  n_labeled_total  <- 0

  for (k in seq_along(conditions)) {
    cnd   <- conditions[[k]]
    cname <- if (!is.null(cnd$name) && nzchar(cnd$name)) as.character(cnd$name) else paste0("Condition", k)
    ccol  <- cnd$col
    cmods <- cnd$modalities

    if (is.null(ccol) || length(ccol) != 1 || !nzchar(ccol)) {
      stop(sprintf("Condition '%s' has no column selected.", cname))
    }
    if (!ccol %in% colnames(clinic_annot)) {
      stop(sprintf("Condition column '%s' not found in clinic_annot.", ccol))
    }

    cvals <- trim_chr(clinic_annot[[ccol]])
    cids  <- ids[cvals %in% trim_chr(cmods)]
    if (length(cids) == 0) {
      stop(sprintf("Condition '%s' matched no sample.", cname))
    }

    meta_k <- data.frame(response = map_response(cids), row.names = cids,
                         stringsAsFactors = FALSE)
    n_labeled_total <- n_labeled_total + sum(!is.na(meta_k$response))

    cond_names            <- c(cond_names, cname)
    cond_scores_list[[k]] <- scores_df[cids, , drop = FALSE]
    cond_meta_list[[k]]   <- meta_k
    cond_ids_list[[k]]    <- cids
  }

  cond_names <- make.unique(cond_names)  # keep condition labels unique for the plots

  if (n_labeled_total == 0) {
    stop("No sample carries a responder / non-responder label after mapping; check the response column and modalities.")
  }

  # Per-condition clinical-value frames (samples x continuous variables).
  cond_clin_list <- NULL
  if (!is.null(clin_num) && length(cont_vars) > 0) {
    cond_clin_list <- lapply(cond_ids_list, function(cids) clin_num[cids, , drop = FALSE])
  }

  ## ---- Run the Python panels ----
  report(0.7, "building boxplots, correlation matrices and ROC curves")
  res <- run_signatures_comparison(
    cond_scores        = cond_scores_list,
    cond_metas         = cond_meta_list,
    cond_names         = as.list(cond_names),
    response_col       = "response",
    responder_label    = responder_label,
    nonresponder_label = nonresponder_label,
    corr_method        = corr_method,
    corr_fdr           = corr_fdr,
    cond_clin          = cond_clin_list,
    out_dir            = output_path
  )

  report(0.95, "finalizing")

  # Normalize the per-subset correlation entries into a tidy R list.
  corr_list <- lapply(res$corr, function(entry) {
    list(label = as.character(entry$label),
         path  = as.character(entry$path),
         n     = as.integer(entry$n))
  })

  clinical_boxplot_path <- if (is.null(res$clinical_boxplot)) NULL else as.character(res$clinical_boxplot)

  list(
    boxplot_path          = as.character(res$boxplot),
    roc_path              = as.character(res$roc),
    corr                  = corr_list,
    stats                 = as.data.frame(res$stats),
    clinical_boxplot_path = clinical_boxplot_path,
    scores                = scores_df,
    output_path           = output_path
  )
}
