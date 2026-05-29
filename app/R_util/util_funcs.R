### This file contains diverse functions useful 
# doDGE()  (DESeq2)
# gsea_multi()
# cox_on_sig_genes()  (cox model KM plot)


doDGEv2 <- function(rnamat,
                    annot,
                    design,
                    condition,
                    modalities_control,
                    modalities_test) {
  
  library(DESeq2)
  
  
  ## Checks
  if (!condition %in% colnames(annot)) {
    stop(sprintf("condition column '%s' not found in annot.", condition))
  }
  
  all_modalities <- as.character(annot[[condition]])
  g1 <- as.character(modalities_control)
  g2 <- as.character(modalities_test)
  
  if (length(intersect(g1, g2)) > 0) {
    stop("Some modalities are present in both modalities_control and modalities_test.")
  }
  
  if (!all(all_modalities %in% c(g1, g2))) {
    missing_modalities <- unique(all_modalities[!all_modalities %in% c(g1, g2)])
    stop(sprintf(
      "Some modalities in annot[['%s']] are not assigned to either group: %s",
      condition,
      paste(missing_modalities, collapse = ", ")
    ))
  }
  
  ## Collapse condition into 2 groups
  annot[[condition]] <- ifelse(all_modalities %in% g1, "control", "test")
  annot[[condition]] <- factor(annot[[condition]], levels = c("control", "test"))
  
  ## Build DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = rnamat,
    colData = annot,
    design = design
  )
  
  dds <- DESeq(dds)
  res <- results(dds, contrast = c(condition, "control", "test"))
  
  res <- as.data.frame(res)
  res <- res[!is.na(res$stat), ]
  
  unique_res <- res[!duplicated(rownames(res)), ]
  significant_genes <- unique_res$stat
  names(significant_genes) <- unique_res$gs
  
  return(dds)
}



# DESeq2
doDGE <- function(rnamat, annot, design){
  library(DESeq2)

  if (!is.null(grouping)) {
    annot$condition <- unname(grouping[as.character(annot$condition)])
    annot$condition <- factor(annot$condition)
  }

    dds <- DESeqDataSetFromMatrix(countData = rnamat, 
                                  colData = annot, 
                                  design = design)  ### attention dans les noms de colonne de annot il faut la colonne "condition"

	
    dds <- DESeq(dds)
    res <- results(dds)

	res = as.data.frame(res)
	res = res[!is.na(res$stat), ]

	unique_res = res[!duplicated(rownames(res)), ]
	significant_genes <- unique_res$stat
	names(significant_genes) = unique_res$gs

    return(dds)
}

# GSEA
gsea_multi <- function(vec_gene, pathw_to_use, min_size = 1) {
  library(fgsea)

  res <- fgseaMultilevel(pathw_to_use, vec_gene, minSize = min_size, BPPARAM = BiocParallel::MulticoreParam(6))#, scoreType = "pos")
  lE_list = res$leadingEdge
  lE_vector = sapply(lE_list, paste, collapse=", ")
  res$leadingEdge = lE_vector
  res = as.data.frame(na.omit(res))
  print(head(res))
  return(res)
}


compute_signature_score <- function(expr_mat, weight_list,
                                    na_as_zero = TRUE,
                                    scale_by_sum_abs = FALSE,
                                    return_df = TRUE,
                                    id_col = "ID_Patient") {
  # expr_mat: matrix/data.frame genes x samples (rownames = genes, colnames = samples)
  # weight_list: named list or named vector gene -> weight
  
  if (is.data.frame(expr_mat)) expr_mat <- as.matrix(expr_mat)
  stopifnot(!is.null(rownames(expr_mat)))
  stopifnot(!is.null(colnames(expr_mat)))

  # Convert weights to named numeric vector
  weights <- unlist(weight_list, use.names = TRUE)
  weights <- as.numeric(weights)
  names(weights) <- names(unlist(weight_list, use.names = TRUE))

  common_genes <- intersect(names(weights), rownames(expr_mat))
  if (length(common_genes) == 0) {
    stop("No overlapping genes between weights and expression matrix.")
  }

  w <- weights[common_genes]
  X <- expr_mat[common_genes, , drop = FALSE]

  # Make sure numeric
  storage.mode(X) <- "numeric"

  if (na_as_zero) {
    X[is.na(X)] <- 0
  }

  scores <- as.numeric(t(w) %*% X)
  names(scores) <- colnames(expr_mat)

  # Optional normalization by total absolute weight
  if (scale_by_sum_abs) {
    denom <- sum(abs(w))
    if (denom > 0) scores <- scores / denom
  }

  if (!return_df) {
    return(scores)  # named numeric vector
  }

  out <- data.frame(
    score = unname(scores),
    stringsAsFactors = FALSE,
    row.names = names(scores)
  )
  out[[id_col]] <- rownames(out)
  out <- out[, c(id_col, "score"), drop = FALSE]

  return(out)
}

Kaplan_Meier_plot <- function(clinic_annot, 
                         subgroup_by, 
                         output_dir,
                         scores_df = NULL,
                         rnaseq = NULL, 
                         survival_time_col = "delpfs",
                         event_realization_col = "PFS",
                         group_quantile = "median", 
                         sample_ID_col = NULL) {
  library(survival)
  library(glmnet)
  library(survminer)


  # Load data
  if (is.null(scores_df)){
    if (!subgroup_by %in% colnames(clinic_annot)){
      if (!subgroup_by %in% rownames(rnaseq)){
        stop(sprintf("Column '%s' not found in clinic_annot and rnaseq genes", subgroup_by))
      }else{
        rnaseq = rnaseq[, intersect(colnames(rnaseq), rownames(clinic_annot))]
        clinic_annot = clinic_annot[colnames(rnaseq), ]
        clinic_annot[[subgroup_by]] <- as.numeric(unlist(rnaseq[subgroup_by, ]))
      }
    }
  }else{
    stopifnot(is.data.frame(scores_df), "score" %in% names(scores_df))
    sc <- as.numeric(scores_df$score)
    names(sc) <- rownames(scores_df)

    # map scores onto clinic_annot by sample_ID_col
    stopifnot(!is.null(sample_ID_col), sample_ID_col %in% names(clinic_annot))
    clinic_annot[[subgroup_by]] <- sc[ as.character(clinic_annot[[sample_ID_col]]) ]
  }

  # Filter out NA on selected columns
  clinicannot_noNA <- clinic_annot[
    !is.na(clinic_annot[[survival_time_col]]) &
    !is.na(clinic_annot[[event_realization_col]]) &
    !is.na(clinic_annot[[subgroup_by]]) &
    !is.na(clinic_annot[[sample_ID_col]]),
    ,
    drop = FALSE
  ]

  # Create plotting groups depending on variable type
  x <- clinicannot_noNA[[subgroup_by]]

  if (is.numeric(x) || is.integer(x)) {
    
    if (group_quantile == "median") {
      med <- median(x, na.rm = TRUE)
      clinicannot_noNA$KM_group <- ifelse(x <= med, "Lower median", "Upper median")
      clinicannot_noNA$KM_group <- factor(clinicannot_noNA$KM_group,
                                          levels = c("Lower median", "Upper median"))
      
    } else if (group_quantile == "tertile") {
      qs <- quantile(x, probs = c(1/3, 2/3), na.rm = TRUE, type = 7)
      clinicannot_noNA$KM_group <- cut(
        x,
        breaks = c(-Inf, qs[1], qs[2], Inf),
        labels = c("T1", "T2", "T3"),
        include.lowest = TRUE
      )
      clinicannot_noNA$KM_group <- factor(clinicannot_noNA$KM_group,
                                          levels = c("T1", "T2", "T3"))
      
    } else if (group_quantile == "quartile") {
      qs <- quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, type = 7)
      clinicannot_noNA$KM_group <- cut(
        x,
        breaks = c(-Inf, qs[1], qs[2], qs[3], Inf),
        labels = c("Q1", "Q2", "Q3", "Q4"),
        include.lowest = TRUE
      )
      clinicannot_noNA$KM_group <- factor(clinicannot_noNA$KM_group,
                                          levels = c("Q1", "Q2", "Q3", "Q4"))
      
    } else {
      stop("group_quantile must be one of: 'median', 'tertile', 'quartile'")
    }
    
  } else {
    clinicannot_noNA$KM_group <- factor(as.character(x))
  }

  # Drop empty levels if any
  clinicannot_noNA$KM_group <- droplevels(clinicannot_noNA$KM_group)

  # Build survival model
  # fit <- survfit(
  #   as.formula(
  #     paste0("Surv(", survival_time_col, ", ", event_realization_col, ") ~ KM_group")
  #   ),
  #   data = clinicannot_noNA
  # )

  clinicannot_noNA$.__time  <- clinicannot_noNA[[survival_time_col]]
  clinicannot_noNA$.__event <- clinicannot_noNA[[event_realization_col]]

  fit <- survfit(
    Surv(.__time, .__event) ~ KM_group,
    data = clinicannot_noNA
  )

  # Plot
  p <- ggsurvplot(
    fit,
    data = clinicannot_noNA,
    pval = TRUE,
    conf.int = TRUE,
    conf.int.style = "step",
    risk.table = "absolute",
    risk.table.height = 0.25,
    risk.table.title = "Number at risk",
    xlab = "Time",
    ylab = "Survival probability",
    legend.title = subgroup_by,
    legend.labs = levels(clinicannot_noNA$KM_group)
  )

  # Save full figure (plot + risk table)
  save_filepath <- file.path(output_dir, paste0("KM_plot_", subgroup_by, "_", group_quantile, ".png"))
  png(
    filename = save_filepath,
    width = 2500,
    height = 2500,
    res = 500
  )
  print(p)
  dev.off()
  return(list(plot = p, save_path = save_filepath))
}

### A décommenter plus tard pour intégrer tout ça au shiny ###

# # Kaplan Meier + fitting of hazard ratios etc....
# cox_on_sig_genes <- function(expr, sig_genes, clin,
#                              time_col, status_col, id_col = NULL,
#                              additional_metric = NULL,
#                              penalty = c("none","lasso","ridge"),
#                              group_rule = c("median","tertile","quartile"),
#                              alpha = 1,        # only used if penalty == "lasso" or custom elastic net
#                              seed = 1, verbose = TRUE) {
#   library(survival)
#   library(glmnet)
#   library(survminer)

#   stopifnot(time_col %in% names(clin), status_col %in% names(clin))
#   penalty   <- match.arg(penalty)
#   group_rule <- match.arg(group_rule)

#   # 0) Coerce expr to matrix, ensure rownames are genes
#   expr <- as.matrix(expr)
#   rn <- rownames(expr); cn <- colnames(expr)
#   if (is.null(rn) || is.null(cn)) stop("expr must have rownames=genes and colnames=samples.")
#   # Try to detect orientation: if more sig_genes match columns than rows, transpose.
#   if (sum(sig_genes %in% rn) < sum(sig_genes %in% cn)) {
#     if (verbose) message("Transposing expr to make rows=genes, cols=samples.")
#     expr <- t(expr)
#     rn <- rownames(expr); cn <- colnames(expr)
#   }

#   # 1) Match samples between expr and clin
#   if (is.null(id_col)) {
#     # try to infer: intersect of colnames(expr) with a column in clin
#     inter_candidates <- intersect(cn, unlist(clin[ , sapply(clin, is.character) | sapply(clin, is.factor), drop=FALSE]))
#     if (length(inter_candidates) > 0) {
#       stop("Please provide id_col (no unique sample ID column was inferred).")
#     } else {
#       stop("id_col is required to align samples (no inference possible).")
#     }
#   }
#   stopifnot(id_col %in% names(clin))
#   clin_ids <- as.character(clin[[id_col]])
#   keep <- intersect(cn, clin_ids)
#   if (length(keep) < 5) stop("Few or no overlapping samples between expression and clinical table.")
#   expr <- expr[, keep, drop = FALSE]
#   clin <- clin[match(keep, clin_ids), , drop = FALSE]
#   rownames(clin) <- clin[[id_col]]

#   # 2) Clean survival columns
#   time  <- suppressWarnings(as.numeric(clin[[time_col]]))
#   status_raw <- clin[[status_col]]
#   if (is.logical(status_raw)) {
#     status <- as.integer(status_raw)  # TRUE=1 event
#   } else if (is.numeric(status_raw)) {
#     status <- as.integer(status_raw)
#   } else {
#     status_map <- c("dead"=1,"deceased"=1,"event"=1,"1"=1,
#                     "alive"=0,"censored"=0,"0"=0)
#     status <- tolower(as.character(status_raw))
#     status <- unname(ifelse(status %in% names(status_map), status_map[status], NA))
#     status <- as.integer(status)
#   }
#   if (anyNA(time) || anyNA(status)) stop("time/status contain NA after coercion. Clean your clinical columns.")
#   surv_obj <- survival::Surv(time, status)

#   # 3) Subset to significant genes present in expr
#   present <- intersect(sig_genes, rownames(expr))
#   if (length(present) < 2) stop("Fewer than 2 significant genes found in the expression matrix.")
#   X <- t(expr[present, , drop = FALSE])      # samples x genes
#   # Drop zero-variance predictors
#   nzv <- apply(X, 2, function(z) sd(z, na.rm=TRUE) > 0)
#   X <- X[, nzv, drop = FALSE]
#   present <- colnames(X)
#   if (length(present) < 2) stop("All selected genes are zero-variance; nothing to model.")
#   # Z-score per gene
#   X <- scale(X)

#   # 4) Fit model
#   set.seed(seed)
#   if (penalty == "none") {
#     # If p is large vs n, this may overfit. Consider lasso/ridge.
#     df <- data.frame(time = time, status = status, X, check.names = FALSE)
#     fml <- as.formula(paste0("survival::Surv(time, status) ~ ", paste(colnames(X), collapse = " + ")))
#     fit <- survival::coxph(fml, data = df, ties = "efron", x = TRUE, y = TRUE)
#     risk <- as.numeric(stats::predict(fit, type = "lp"))
#     coefs <- stats::coef(fit)
#     # PH test
#     zph  <- try(survival::cox.zph(fit), silent = TRUE)
#   } else {
#     # Penalized Cox via glmnet
#     if (!requireNamespace("glmnet", quietly = TRUE))
#       stop("Package 'glmnet' is required for penalized Cox.")
#     a <- if (penalty == "lasso") 1 else 0 # ridge=0
#     if (!missing(alpha) && penalty != "lasso" && penalty != "ridge") a <- alpha
#     cvfit <- glmnet::cv.glmnet(X, surv_obj, family = "cox", alpha = a, nfolds = 10)
#     fit   <- cvfit
#     risk  <- as.numeric(predict(cvfit, newx = X, s = "lambda.min", type = "link"))
#     # Extract nonzero coefficients at lambda.min
#     b <- as.matrix(glmnet::coef.glmnet(cvfit$glmnet.fit, s = cvfit$lambda.min))
#     nz <- which(abs(b) > 0); coefs <- b[nz, , drop = TRUE]; names(coefs) <- rownames(b)[nz]
#     # Optionally refit an unpenalized Cox on selected features (for zph, HR/CI):
#     if (length(coefs) > 0 && length(coefs) <= min(30, nrow(X)-5)) {
#       # Selected features
#       vars <- names(coefs)
#       stopifnot(length(vars) > 0)

#       # Keep original (possibly non-syntactic) names
#       df2 <- data.frame(time = time, status = status,
#                         X[, vars, drop = FALSE],
#                         check.names = FALSE)

#       # Backtick every term so symbols like "-" or ":" are safe
#       bt   <- sprintf("`%s`", vars)
#       fml2 <- as.formula(paste0("survival::Surv(time, status) ~ ", paste(bt, collapse = " + ")))

#       refit <- survival::coxph(fml2, data = df2, ties = "efron", x = TRUE, y = TRUE)
#       zph   <- try(survival::cox.zph(refit), silent = TRUE)
#     } else {
#       zph <- NULL
#     }
#   }

#   names(risk) <- rownames(X)

#   # 5) Stratify and KM
#   cut <- switch(group_rule,
#                 median   = stats::median(risk),
#                 tertile  = as.numeric(stats::quantile(risk, probs = 2/3, na.rm=TRUE)),
#                 quartile = as.numeric(stats::quantile(risk, probs = 0.75, na.rm=TRUE)))
#   group <- if (group_rule == "median") ifelse(risk > cut, "High", "Low") else
#            ifelse(risk >= cut, "High", "Low")
#   if (is.null(additional_metric)) {
#      group <- factor(group, levels = c("Low","High"))
#   }
#   else{
#     group <- factor(paste0(group,clin[,additional_metric]))
#     print(levels(group))
#   }
  
#   df_plot <- data.frame(time = time, status = status, group = group)
#   km <- survival::survfit(survival::Surv(time, status) ~ group, data = df_plot)

#   # 6) Concordance (C-index)
#   cidx <- try({
#     survival::concordance(surv_obj ~ risk)$concordance
#   }, silent = TRUE)
#   if (inherits(cidx, "try-error")) cidx <- NA_real_

#   # 7) Plot KM (uses survminer if available)
#   plot_obj <- NULL
#   if (requireNamespace("survminer", quietly = TRUE)) {
#     df_plot <- data.frame(time = time, status = status, group = group)
#     plot_obj <- survminer::ggsurvplot(km, data = df_plot,
#                                       risk.table = TRUE, conf.int = TRUE,
#                                       pval = TRUE, legend.title = "Risk",
#                                       legend.labs = levels(group))
#   } else {
#     if (verbose) {
#       plot(km, col = c("black","red"), lwd = 2, mark.time = TRUE,
#            xlab = "Time", ylab = "Survival probability", main = "KM by risk group")
#       legend("bottomleft", legend = levels(group), col = c("black","red"), lwd = 2)
#     }
#   }

#   if (verbose) {
#     message("# Samples: ", nrow(X),
#             " | # Genes used: ", ncol(X),
#             " | Penalty: ", penalty,
#             " | C-index: ", round(cidx, 3))
#   }

#   invisible(list(
#     penalty    = penalty,
#     genes_used = colnames(X),
#     model      = fit,
#     coefs      = coefs,
#     risk_score = risk,
#     cut_value  = cut,
#     group_rule = group_rule,
#     groups     = group,
#     km         = km,
#     c_index    = cidx,
#     zph        = zph,
#     plot       = plot_obj
#   ))
# }


