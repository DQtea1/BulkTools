# R_analysis/ICA.R
# R wrappers around py/py_ICA.py for the ICA module (mod_ICA.R).
#   mstd_pipe() : part 1 - Most Stable Transcriptomic Dimension (sica.base.MSTD).
#   ica_pipe()  : part 2 - stabilized ICA, saves A/S (CSV + RDS), builds all plots.
# Follows the user's notebooks (04_PCA_ICA/02_ICA.ipynb): X is samples x genes,
#   A = sICA.S_ (metagenes x genes), S = pinv(A).T @ X.T (metagenes x samples).
# Python is driven through reticulate, like R_analysis/signature_projection.R.

# Coerce to a numeric matrix (genes x samples) and optionally VST-normalize.
.ica_prepare_expr <- function(rnafilt_counts, do_vst) {
    expr <- as.matrix(rnafilt_counts)
    if (isTRUE(do_vst)) {
        expr <- normVST_bulk(expr)
    }
    storage.mode(expr) <- "double"
    expr
}

mstd_pipe <- function(rnafilt_counts, output_dir, m, M, step, n_runs,
                      do_vst = TRUE, progress_cb = NULL) {
    library(reticulate)

    report <- function(frac, detail) {
        if (is.function(progress_cb)) progress_cb(frac, detail)
    }

    use_python("/opt/conda/envs/BulkTools/bin/python", required = TRUE)
    source_python("py/py_ICA.py")

    report(0.1, "preprocessing (optional VST)")
    expr <- .ica_prepare_expr(rnafilt_counts, do_vst)

    out_dir <- file.path(output_dir, "ICA")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    report(0.4, "running MSTD over the dimension range (this can take a while)")
    res <- compute_mstd(
        expr         = expr,
        gene_names   = rownames(expr),
        sample_names = colnames(expr),
        m            = as.integer(m),
        M            = as.integer(M),
        step         = as.integer(step),
        n_runs       = as.integer(n_runs),
        out_dir      = out_dir
    )

    report(1.0, "done")
    list(mstd_plot = res$save_path)
}

ica_pipe <- function(rnafilt_counts, output_dir, n_components, n_runs,
                     do_vst = TRUE, clinic_annot = NULL,
                     scatter_x = NULL, scatter_y = NULL, color_col = NULL,
                     scatter_bins = 15, progress_cb = NULL) {
    library(reticulate)

    report <- function(frac, detail) {
        if (is.function(progress_cb)) progress_cb(frac, detail)
    }

    use_python("/opt/conda/envs/BulkTools/bin/python", required = TRUE)
    source_python("py/py_ICA.py")

    report(0.1, "preprocessing (optional VST)")
    expr <- .ica_prepare_expr(rnafilt_counts, do_vst)
    gene_names   <- rownames(expr)
    sample_names <- colnames(expr)

    out_dir <- file.path(output_dir, "ICA")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    report(0.35, "fitting stabilized ICA")
    fit <- run_stabilized_ica(
        expr         = expr,
        gene_names   = gene_names,
        sample_names = sample_names,
        n_components = as.integer(n_components),
        n_runs       = as.integer(n_runs)
    )

    comp_names <- unlist(fit$comp_names)
    K <- length(comp_names)

    # A = metagenes x genes (gene weights) ; S = metagenes x samples (activities).
    A_df <- as.data.frame(fit$A)
    rownames(A_df) <- comp_names
    colnames(A_df) <- gene_names

    S_df <- as.data.frame(fit$S)
    rownames(S_df) <- comp_names
    colnames(S_df) <- sample_names

    stability_df <- data.frame(
        component       = comp_names,
        stability_index = as.numeric(fit$stability),
        row.names       = comp_names,
        check.names     = FALSE
    )

    report(0.6, "saving A / S matrices (CSV + RDS)")
    a_csv <- file.path(out_dir, sprintf("A_matrix_%dmetagenes.csv", K))
    s_csv <- file.path(out_dir, sprintf("S_matrix_%dmetagenes.csv", K))
    write_csv_mkdir(A_df, a_csv)
    write_csv_mkdir(S_df, s_csv)
    write_csv_mkdir(stability_df, file.path(out_dir, "stability_index.csv"), row.names = FALSE)
    saveRDS(A_df, sub("\\.csv$", ".rds", a_csv))
    saveRDS(S_df, sub("\\.csv$", ".rds", s_csv))

    # Default scatter metagenes: first two if not provided / invalid.
    if (is.null(scatter_x) || !scatter_x %in% comp_names) scatter_x <- comp_names[1]
    if (is.null(scatter_y) || !scatter_y %in% comp_names) {
        scatter_y <- if (K >= 2) comp_names[2] else comp_names[1]
    }
    if (!is.null(color_col) && !nzchar(color_col)) color_col <- NULL

    report(0.8, "building plots")
    scatter_path <- plot_scatter_marginals(
        S = fit$S, comp_names = comp_names, sample_names = sample_names,
        comp_x = scatter_x, comp_y = scatter_y, out_dir = out_dir,
        clinic = clinic_annot, color_col = color_col, bins = as.integer(scatter_bins)
    )
    activity_heatmap_path <- plot_activity_heatmap(fit$S, comp_names, sample_names, out_dir)
    source_dist_path      <- plot_source_distributions(fit$A, comp_names, gene_names, out_dir)
    stability_path        <- plot_stability_index(fit$stability, comp_names, out_dir)
    corr_path             <- plot_component_correlation(fit$S, comp_names, out_dir)

    report(0.92, "clinical association (if clinic provided)")
    clinical <- plot_clinical_association(fit$S, comp_names, sample_names, clinic_annot, out_dir)

    report(1.0, "done")
    list(
        scatter          = scatter_path,
        activity_heatmap = activity_heatmap_path,
        source_dist      = source_dist_path,
        stability        = stability_path,
        corr             = corr_path,
        clinical         = clinical,
        A                = A_df,
        S                = S_df,
        comp_names       = comp_names,
        out_dir          = out_dir
    )
}
