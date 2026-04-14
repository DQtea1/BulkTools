

signature_proj_pipe = function(rnafilt_counts, clinic_annot, 
                                output_dir, therapy_used, signature_name,
                                signature_to_use, sample_to_project_path, 
                                contrast, resp_var, non_resp_var, clinic_filters = NULL,
                                filter_by_gene, keep_low_or_high, quantile_thr,
                                survival_time_col = "delpfs", event_realization_col = "PFS", 
                                sample_ID_col = "ID_Patient", group_quantile = "median"){
    library(reticulate)

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

    prepared_data = prepare_projection(rnafilt_counts, 
                                       sample_to_project_path)

    merged_bulks_vst = prepared_data$merged_bulks_vst
    ID_ref_samples   = prepared_data$ID_ref_samples
    sample_name      = prepared_data$sample_name

    reference_scores = compute_signature_score(merged_bulks_vst[, ID_ref_samples, drop = FALSE], signature_to_use)

    sample_score = compute_signature_score(merged_bulks_vst[, sample_name, drop = FALSE], signature_to_use)
    
    # 1st panel : Quantile + sample projection + accuracy metrics in python :
    clean_ids = function(x) {
    x = trimws(as.character(x))
    x = x[!is.na(x)]
    x = x[x != ""]
    x = x[x != "NA"]
    unique(x)
    }

    responder_ids     = clean_ids(clinic_annot[clinic_annot[[contrast]] == resp_var, "ID_Patient"])
    non_responder_ids = clean_ids(clinic_annot[clinic_annot[[contrast]] == non_resp_var, "ID_Patient"])
    responder_ids     = intersect(responder_ids, ID_ref_samples)
    non_responder_ids = intersect(non_responder_ids, ID_ref_samples)

    model_eval_dir = file.path(output_dir, "model_eval")
    dir.create(model_eval_dir, recursive = TRUE, showWarnings = FALSE)


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

    # 3rd panel Signature evaluation (ROC curve, boxplot, conf mat? ) in python :
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
