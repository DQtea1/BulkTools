

signature_proj_pipe = function(rnafilt_counts, clinic_annot,
                                output_dir, therapy_used, signature_name,
                                signature_to_use, sample_to_project_path,
                                contrast, resp_var, non_resp_var, clinic_filters = NULL,
                                filter_by_gene, keep_low_or_high, quantile_thr,
                                survival_time_col = "delpfs", event_realization_col = "PFS",
                                sample_ID_col = "ID_Patient", group_quantile = "median",
                                do_km_plot = TRUE, extra_components = NULL,
                                progress_cb = NULL){
    library(reticulate)

    # Optional progress reporter: progress_cb(frac, detail), frac in [0, 1].
    report <- function(frac, detail) {
        if (is.function(progress_cb)) progress_cb(frac, detail)
    }
    report(0.1, "preprocessing & clinic/gene filtering")

    # Combine the base signature with optional extra signatures / single genes
    # at the SCORE level: each component's per-sample score is computed, then
    # added / subtracted / multiplied / divided into the running score.
    compute_combined_score <- function(expr_mat, base_geneset, components = NULL) {
        combined <- compute_signature_score(expr_mat, base_geneset, scale_by_sum_abs = TRUE)
        if (!is.null(components) && length(components) > 0) {
            for (comp in components) {
                sc_i <- compute_signature_score(expr_mat, comp$geneset, scale_by_sum_abs = TRUE)
                sc_i <- sc_i[rownames(combined), , drop = FALSE]
                combined$score <- switch(comp$op,
                    "+" = combined$score + sc_i$score,
                    "-" = combined$score - sc_i$score,
                    "*" = combined$score * sc_i$score,
                    "/" = combined$score / sc_i$score,
                    stop(sprintf("Unsupported combination operator '%s'.", comp$op))
                )
            }
        }
        combined
    }

    use_python("/opt/conda/envs/BulkTools/bin/python", required = TRUE)
    source_python("py/py_plots.py")
    source_python("py/py_models.py")

    parsed_design = "nul"

    # Filtering samples on selected gene's expression and clinic features
    filtered_data = filtering_on_clinic_and_genes(clinic_annot, rnafilt_counts, 
                                                  clinic_filters     = clinic_filters,
                                                  filter_by_gene     = filter_by_gene, 
                                                  keep_low_or_high   = keep_low_or_high, 
                                                  quantile_thr       = quantile_thr,
                                                  parsed_design      = parsed_design,
                                                  output_dir         =  output_dir,
                                                  folder_name        = "signature_proj"
                                                  ) 
    
    rnafilt_counts = filtered_data$rnafilt_counts
    clinic_annot   = filtered_data$clinic_annot
    output_DESeq   = filtered_data$output_path

    report(0.3, "importing projected sample + VST normalization")
    prepared_data = prepare_projection(rnafilt_counts,
                                       sample_to_project_path)

    merged_bulks_vst = prepared_data$merged_bulks_vst
    ID_ref_samples   = prepared_data$ID_ref_samples
    sample_name      = prepared_data$sample_name

    report(0.45, "computing (combined) signature scores")
    reference_scores = compute_combined_score(merged_bulks_vst[, ID_ref_samples, drop = FALSE], signature_to_use, extra_components)

    sample_score = compute_combined_score(merged_bulks_vst[, sample_name, drop = FALSE], signature_to_use, extra_components)
    
    # 1st panel : Quantile + sample projection + accuracy metrics in python :
    clean_ids = function(x) {
    x = trimws(as.character(x))
    x = x[!is.na(x)]
    x = x[x != ""]
    x = x[x != "NA"]
    unique(x)
    }

    if (is.null(contrast) || length(contrast) == 0 || !nzchar(contrast)) {
        stop("Please select a single clinical column under 'Contrast'.")
    }
    if (length(contrast) > 1) {
        stop("Only one clinical column can be used as 'Contrast'.")
    }
    if (!contrast %in% colnames(clinic_annot)) {
        stop(sprintf("Contrast column '%s' not found in clinic_annot.", contrast))
    }
    if (is.null(resp_var) || length(resp_var) == 0 || all(!nzchar(resp_var))) {
        stop("Please select at least one modality under 'Responders'.")
    }
    if (is.null(non_resp_var) || length(non_resp_var) == 0 || all(!nzchar(non_resp_var))) {
        stop("Please select at least one modality under 'Non-Responders'.")
    }

    if (!"ID_Patient" %in% colnames(clinic_annot)) {
        if (!is.null(rownames(clinic_annot))) {
            clinic_annot$ID_Patient <- rownames(clinic_annot)
        } else {
            stop("clinic_annot has no 'ID_Patient' column and no rownames; cannot identify samples.")
        }
    }

    clinic_ids_all    = trimws(as.character(clinic_annot$ID_Patient))
    ref_ids_all       = trimws(as.character(ID_ref_samples))
    contrast_values   = trimws(as.character(clinic_annot[[contrast]]))
    resp_var_trim     = trimws(as.character(resp_var))
    non_resp_var_trim = trimws(as.character(non_resp_var))

    resp_mask         = contrast_values %in% resp_var_trim
    non_resp_mask     = contrast_values %in% non_resp_var_trim

    responder_ids_raw     = clean_ids(clinic_ids_all[resp_mask])
    non_responder_ids_raw = clean_ids(clinic_ids_all[non_resp_mask])
    responder_ids         = intersect(responder_ids_raw,     ref_ids_all)
    non_responder_ids     = intersect(non_responder_ids_raw, ref_ids_all)

    if (length(responder_ids) < 2 || length(non_responder_ids) < 2) {
        unique_modalities <- sort(unique(contrast_values))
        clinic_overlap    <- intersect(clinic_ids_all, ref_ids_all)
        stop(sprintf(
            paste0(
                "Need at least 2 samples per class for nested CV.\n",
                "--- DIAGNOSTIC ---\n",
                "clinic_annot rows after filter        : %d\n",
                "contrast column                       : '%s'\n",
                "  unique values (first 15)            : %s\n",
                "  rows where contrast %%in%% [%s]      : %d\n",
                "  rows where contrast %%in%% [%s]      : %d\n",
                "reference bulk sample count           : %d\n",
                "clinic ID_Patient ∩ ref bulk          : %d / %d\n",
                "  clinic IDs (first 5)                : %s\n",
                "  ref bulk IDs (first 5)              : %s\n",
                "responder IDs intersected with bulk   : %d -> [%s]\n",
                "non-responder IDs intersected w/ bulk : %d -> [%s]\n",
                "------------------\n",
                "Common causes: (a) selected modalities do not appear verbatim in the contrast column, ",
                "(b) clinic ID_Patient values don't match the bulk column names (whitespace, separator, case), ",
                "(c) clinic_filter / gene_filter dropped all responder or non-responder samples."
            ),
            nrow(clinic_annot),
            contrast,
            paste(head(unique_modalities, 15), collapse = ", "),
            paste(resp_var_trim, collapse = ", "), sum(resp_mask),
            paste(non_resp_var_trim, collapse = ", "), sum(non_resp_mask),
            length(ref_ids_all),
            length(clinic_overlap), length(clinic_ids_all),
            paste(head(clinic_ids_all, 5), collapse = ", "),
            paste(head(ref_ids_all, 5), collapse = ", "),
            length(responder_ids),     paste(head(responder_ids, 5), collapse = ", "),
            length(non_responder_ids), paste(head(non_responder_ids, 5), collapse = ", ")
        ))
    }

    model_eval_dir = file.path(output_dir, "model_eval")
    dir.create(model_eval_dir, recursive = TRUE, showWarnings = FALSE)


    report(0.55, "nested cross-validation (this can take a while)")
    res_nested_cv = nested_cv_signature(
                                        df_scores                = reference_scores[rownames(reference_scores) %in% ID_ref_samples, ],                         # index = sample IDs
                                        sample_ID_responders     = responder_ids,
                                        sample_ID_non_responders = non_responder_ids,
                                        score_col                = "score",
                                        n_outer                  = 20,
                                        n_inner                  = 20,
                                        threshold_criterion      = "youden",
                                        use_gray_zone            = TRUE,
                                        gray_target_sensitivity  = 0.90,
                                        gray_target_specificity  = 0.90,
                                        n_bootstrap              = 1000,
                                        ci_level                 = 0.95
                                        )

    # Query = projected sample(s), not reference scores
    query_df     = sample_score[, , drop = FALSE]
    query_scores = as.numeric(query_df[["score"]])
    query_ids    = rownames(query_df)

    stopifnot(length(query_scores) == length(query_ids))
    query_ids_py = as.list(unname(query_ids))

    conf_report = sample_confidence_report(
                                           scores                  = unname(query_scores),
                                           res_nested_cv_signature = res_nested_cv,
                                           sample_ids              = query_ids_py,
                                           n_bootstrap             = 1000
                                           )

    sample_proj_dir = file.path(output_dir, "projection")
    dir.create(sample_proj_dir, recursive = TRUE, showWarnings = FALSE)
    save_path = file.path(sample_proj_dir, paste0("proj_", paste(query_ids, collapse = "_"), ".png"))

    report(0.8, "building projection + confidence plot")
    proj_plot = plot_sample_signature_confidence(
                                                 query_scores       = unname(query_scores),
                                                 res_v2             = res_nested_cv,
                                                 query_ids          = query_ids_py,
                                                 confidence_res     = conf_report,
                                                 confidence_kwargs  = list(n_bootstrap = 1000, use_gray_zone = TRUE),
                                                 n_bootstrap_threshold_plot = 1000,
                                                 language           = "en",
                                                 save_path          = save_path,
                                                 figsize            = c(13, 8)
                                                 )


    # 2nd panel Kaplan-Meier plot in R :
    KM_plot = NULL
    if (isTRUE(do_km_plot)) {
        valid_col <- function(x) !is.null(x) && length(x) == 1 && nzchar(x)
        if (!valid_col(survival_time_col) || !survival_time_col %in% colnames(clinic_annot)) {
            stop(sprintf(
                "KM plot requested but 'Survival time column' is missing or not in clinic_annot (got: '%s'). Pick a valid column or set 'Plot KM plot ?' to 'Nein'.",
                if (is.null(survival_time_col) || length(survival_time_col) == 0) "NULL" else as.character(survival_time_col)
            ))
        }
        if (!valid_col(event_realization_col) || !event_realization_col %in% colnames(clinic_annot)) {
            stop(sprintf(
                "KM plot requested but 'Event realization column' is missing or not in clinic_annot (got: '%s'). Pick a valid column or set 'Plot KM plot ?' to 'Nein'.",
                if (is.null(event_realization_col) || length(event_realization_col) == 0) "NULL" else as.character(event_realization_col)
            ))
        }

        KM_plot = Kaplan_Meier_plot(
                                    clinic_annot           = clinic_annot,
                                    subgroup_by            = signature_name,
                                    output_dir             = output_dir,
                                    scores_df              = reference_scores,
                                    survival_time_col      = survival_time_col,
                                    event_realization_col  = event_realization_col,
                                    group_quantile         = group_quantile,
                                    sample_ID_col          = sample_ID_col
                                    )
    }

    # 3rd panel Signature evaluation (ROC curve, boxplot, conf mat? ) in python :
    report(0.9, "survival / ROC / boxplot")
    box_plot = my_Box_Wilcox(
                             df             = reference_scores[rownames(reference_scores) %in% ID_ref_samples, ],   # data.frame avec rownames = IDs
                             responders     = responder_ids,               # vecteur R
                             non_responders = non_responder_ids,
                             out_dir        = model_eval_dir,
                             score_col      = "score",
                             return_mode    = "path"
                             )

    roc_plot = myROC_AUC(
                         df             = reference_scores[rownames(reference_scores) %in% ID_ref_samples, ],
                         responders     = responder_ids,               # vecteur R
                         non_responders = non_responder_ids,
                         out_dir        = model_eval_dir,
                         score_col      = "score",
                        )


    return(list(proj_plot = proj_plot, 
                KM_plot   = KM_plot, 
                roc_plot  = roc_plot, 
                box_plot  = box_plot))
}
