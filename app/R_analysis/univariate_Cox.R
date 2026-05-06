# Wraps a univariate Cox proportional-hazards screen:
# scale features -> per-feature coxph(Surv(time, event) ~ x) -> table + volcano.
# Input rnafilt_counts is features x samples (same convention as the other pipes).

univariate_cox_pipe <- function(rnafilt_counts,
                                clinic_annot,
                                output_dir,
                                survival_time_col = "OS_time",
                                event_realization_col = "DECES",
                                clinic_filters = NULL,
                                filter_by_gene = NULL,
                                keep_low_or_high = NULL,
                                quantile_thr = NULL,
                                top_n_label = 50,
                                label_filter_pattern = NULL,
                                pathways_to_use = "REACTOME_pathways",
                                gsea_pval_thr = 0.01,
                                gsea_top_n = 70) {
  library(survival)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(fgsea)

  parsed_design <- paste(survival_time_col, event_realization_col, sep = "__")

  filtered_data <- filtering_on_clinic_and_genes(
    clinic_annot, rnafilt_counts,
    clinic_filters   = clinic_filters,
    filter_by_gene   = filter_by_gene,
    keep_low_or_high = keep_low_or_high,
    quantile_thr     = quantile_thr,
    parsed_design    = parsed_design,
    output_dir       = output_dir,
    folder_name      = "univariate_Cox"
  )

  rnafilt_counts <- filtered_data$rnafilt_counts
  clinic_annot   <- ensure_clinic_sample_id_col(filtered_data$clinic_annot)
  output_path    <- filtered_data$output_path

  if (!survival_time_col %in% colnames(clinic_annot)) {
    stop(sprintf("Survival time column '%s' not found in clinic_annot.", survival_time_col))
  }
  if (!event_realization_col %in% colnames(clinic_annot)) {
    stop(sprintf("Event column '%s' not found in clinic_annot.", event_realization_col))
  }

  clinic_annot[[survival_time_col]] <- suppressWarnings(
    as.numeric(as.character(clinic_annot[[survival_time_col]]))
  )
  clinic_annot[[event_realization_col]] <- suppressWarnings(
    as.numeric(as.character(clinic_annot[[event_realization_col]]))
  )

  keep <- !is.na(clinic_annot[[survival_time_col]]) &
          !is.na(clinic_annot[[event_realization_col]])
  clinic_annot <- clinic_annot[keep, , drop = FALSE]
  if (nrow(clinic_annot) == 0) {
    stop("No samples remain after dropping rows with missing survival time or event.")
  }

  common <- intersect(clinic_annot$ID_Patient, colnames(rnafilt_counts))
  if (length(common) == 0) {
    stop("No samples shared between rnafilt_counts and clinic_annot after filtering.")
  }
  clinic_annot <- clinic_annot[match(common, clinic_annot$ID_Patient), , drop = FALSE]
  rnafilt_counts <- rnafilt_counts[, common, drop = FALSE]

  surv_obj <- Surv(
    time  = clinic_annot[[survival_time_col]],
    event = clinic_annot[[event_realization_col]]
  )

  # Cox expects samples x features
  expr_mat <- t(as.matrix(rnafilt_counts))
  expr_mat_scaled <- scale(expr_mat)
  feature_ok <- apply(expr_mat_scaled, 2, function(x) {
    !all(is.na(x)) && sd(x, na.rm = TRUE) > 0
  })
  expr_mat_scaled <- expr_mat_scaled[, feature_ok, drop = FALSE]

  if (ncol(expr_mat_scaled) == 0) {
    stop("All features were dropped (constant or NA-only after scaling).")
  }

  cox_one_feature <- function(x) {
    fit <- tryCatch(coxph(surv_obj ~ x), error = function(e) NULL)
    if (is.null(fit)) {
      return(data.frame(
        coef       = NA_real_,
        HR         = NA_real_,
        pval       = NA_real_,
        z          = NA_real_,
        score_test = NA_real_,
        score_pval = NA_real_
      ))
    }
    s <- summary(fit)
    score_stat <- if (!is.null(s$sctest["test"])) s$sctest["test"] else s$sctest["statistic"]
    data.frame(
      coef       = s$coef[1, "coef"],
      HR         = exp(s$coef[1, "coef"]),
      pval       = s$coef[1, "Pr(>|z|)"],
      z          = s$coef[1, "z"],
      score_test = score_stat,
      score_pval = s$sctest["pvalue"]
    )
  }

  res_list <- lapply(seq_len(ncol(expr_mat_scaled)), function(j) {
    cox_one_feature(expr_mat_scaled[, j])
  })
  res <- do.call(rbind, res_list)
  res$gene <- colnames(expr_mat_scaled)
  res <- res %>%
    filter(!is.na(pval)) %>%
    mutate(FDR = p.adjust(pval, method = "BH")) %>%
    arrange(pval)

  if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
  write_csv_mkdir(res, file.path(output_path, "univariate_cox_table.csv"), row.names = FALSE)

  ## VOLCANO PLOT ##
  df <- res %>%
    mutate(
      negLog10P = -log10(as.numeric(pval)),
      color = case_when(
        pval < 0.05 & coef > 0 ~ "High expression = worse",
        pval < 0.05 & coef < 0 ~ "High expression = better",
        TRUE ~ "NS"
      )
    ) %>%
    filter(is.finite(negLog10P))

  n_label <- if (is.null(top_n_label) || is.na(top_n_label)) 0L else max(0L, as.integer(top_n_label))
  top_genes <- df %>% arrange(pval) %>% head(n_label)
  if (!is.null(label_filter_pattern) && nzchar(label_filter_pattern)) {
    top_genes <- top_genes %>% filter(!grepl(label_filter_pattern, gene, perl = TRUE))
  }

  p <- ggplot(df, aes(x = log(HR), y = negLog10P, color = color)) +
    geom_point(alpha = 0.5, size = 2) +
    scale_color_manual(values = c(
      "High expression = better" = "#0d55c2",
      "High expression = worse"  = "#d12727",
      "NS"                       = "grey"
    )) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_text_repel(
      data = top_genes,
      aes(label = gene),
      size = 4,
      show.legend = FALSE,
      max.overlaps = Inf,
      box.padding = 0.5,
      point.padding = 0.5,
      segment.color = "grey50"
    ) +
    theme_minimal() +
    labs(
      title = sprintf("Univariate Cox: Surv(%s, %s) ~ feature",
                      survival_time_col, event_realization_col),
      x = "log(HR)",
      y = "-log10(p-value)"
    ) + 
    xlim(-10, 10)


  ggsave(file.path(output_path, "univariate_cox_volcano.pdf"), p, width = 12, height = 10)
  ggsave(file.path(output_path, "univariate_cox_volcano.png"), p, width = 12, height = 10, dpi = 300)

  ## GSEA on the Cox z-statistic ##
  selected_pathways <- normalize_gene_set_collection(
    get(pathways_to_use),
    fallback_prefix = pathways_to_use
  )

  if (length(selected_pathways) == 0) {
    stop(sprintf("The selected pathway collection '%s' does not contain any usable gene set.", pathways_to_use))
  }

  ranked_genes <- res$z
  names(ranked_genes) <- res$gene
  ranked_genes <- ranked_genes[!is.na(ranked_genes)]
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)

  res_gsea <- fgseaMultilevel(
    pathways = selected_pathways,
    stats    = ranked_genes,
    minSize  = 10,
    maxSize  = 500
  )

  res_gsea <- res_gsea[order(abs(res_gsea$NES), decreasing = TRUE), ]
  res_gsea_flat <- flatten_list_cols(res_gsea)
  write_csv_mkdir(res_gsea_flat, file.path(output_path, "gsea_results.csv"), row.names = FALSE)

  top_n <- max(0L, as.integer(gsea_top_n))
  res_sig <- res_gsea %>% filter(pval < gsea_pval_thr)

  top_pos <- res_sig %>% filter(NES > 0) %>% arrange(desc(NES)) %>% head(top_n)
  top_neg <- res_sig %>% filter(NES < 0) %>% arrange(NES)        %>% head(top_n)
  df_top  <- bind_rows(top_pos, top_neg) %>%
    mutate(
      neglog10p     = -log10(pval),
      size_set      = size,
      label_pathway = sub("^[A-Z0-9]+_", "", pathway)
    )

  if (nrow(df_top) > 0) {
    df_top$label_pathway <- factor(
      df_top$label_pathway,
      levels = df_top$label_pathway[order(df_top$NES)]
    )
  }

  p_gsea <- ggplot(df_top, aes(x = NES, y = label_pathway)) +
    geom_point(aes(size = size_set, color = neglog10p), alpha = 0.9) +
    scale_color_viridis_c(option = "viridis") +
    scale_size(range = c(3, 10)) +
    coord_cartesian(xlim = c(-4, 4)) +
    labs(
      x     = "Normalized Enrichment Score",
      y     = "",
      color = "-log10(pval)",
      size  = "Gene set size",
      title = sprintf("GSEA on Cox z-statistic (%s)", pathways_to_use)
    ) +
    theme_bw() +
    theme(
      axis.text.y     = element_text(size = 5),
      panel.grid.minor = element_blank()
    )

  ggsave(file.path(output_path, "gsea_dotplot.pdf"), p_gsea, width = 12, height = 10)
  ggsave(file.path(output_path, "gsea_dotplot.png"), p_gsea, width = 12, height = 10, dpi = 300)

  return(list(
    plots  = list(volcano = p, gsea_dotplot = p_gsea),
    tables = list(cox_results = res, gsea_results = res_gsea_flat)
  ))
}
