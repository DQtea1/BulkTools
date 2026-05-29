format_clinic_filter_description <- function(clinic_filters) {
    normalized_filters <- normalize_clinic_filters(clinic_filters)

    if (is.null(normalized_filters)) {
        return("none")
    }

    paste(vapply(names(normalized_filters), function(filter_name) {
        values <- paste(sort(unique(normalized_filters[[filter_name]])), collapse = ", ")
        sprintf("%s in {%s}", filter_name, values)
    }, character(1)), collapse = "; ")
}


format_group_gene_filter_description <- function(group_gene_filter) {
    normalized_filters <- normalize_group_gene_filters(group_gene_filter)

    if (is.null(normalized_filters)) {
        return("none")
    }

    paste(vapply(normalized_filters, function(filter_def) {
        sprintf(
            "%s %s q=%s",
            filter_def$gene,
            filter_def$keep_low_or_high,
            format(filter_def$quantile_thr, scientific = FALSE, trim = TRUE)
        )
    }, character(1)), collapse = " AND ")
}


group_gene_filter_note <- function(rnafilt_counts, group_gene_filter) {
    if (is.data.frame(rnafilt_counts)) {
        rnafilt_counts <- as.matrix(rnafilt_counts)
    }

    if (!group_gene_filter$gene %in% rownames(rnafilt_counts)) {
        return(sprintf("gene '%s' is missing from the expression matrix", group_gene_filter$gene))
    }

    if (isTRUE(all.equal(group_gene_filter$quantile_thr, 0)) && identical(group_gene_filter$keep_low_or_high, "low")) {
        return("q=0 with the strict 'low' rule is always empty")
    }

    if (isTRUE(all.equal(group_gene_filter$quantile_thr, 1)) && identical(group_gene_filter$keep_low_or_high, "high")) {
        return("q=1 with the strict 'high' rule is always empty")
    }

    gene_expr <- as.numeric(rnafilt_counts[group_gene_filter$gene, , drop = TRUE])
    gene_expr <- gene_expr[!is.na(gene_expr)]

    if (length(gene_expr) > 0 && length(unique(gene_expr)) == 1) {
        return("expression is constant across the remaining samples")
    }

    NULL
}


diagnose_group_matching <- function(clinic_annot,
                                    rnafilt_counts,
                                    clinic_filters = NULL,
                                    group_gene_filter = NULL,
                                    sample_id_col = "ID_Patient") {
    normalized_clinic_filters <- normalize_clinic_filters(clinic_filters)
    normalized_gene_filters <- normalize_group_gene_filters(group_gene_filter)
    available_ids <- intersect(
        unique(extract_clinic_sample_ids(clinic_annot, sample_id_col = sample_id_col)),
        colnames(rnafilt_counts)
    )

    clinic_ids <- available_ids
    if (!is.null(normalized_clinic_filters)) {
        clinic_ids <- intersect(
            available_ids,
            matching_clinic_sample_ids(
                clinic_annot = clinic_annot,
                clinic_filters = normalized_clinic_filters,
                sample_id_col = sample_id_col
            )
        )
    }

    gene_ids <- available_ids
    gene_filter_counts <- list()
    if (!is.null(normalized_gene_filters)) {
        for (filter_idx in seq_along(normalized_gene_filters)) {
            gene_filter <- normalized_gene_filters[[filter_idx]]
            matched_gene_ids <- intersect(
                available_ids,
                matching_gene_sample_ids(rnafilt_counts, gene_filter)
            )

            gene_ids <- intersect(gene_ids, matched_gene_ids)
            gene_filter_counts[[filter_idx]] <- list(
                filter = gene_filter,
                match_count = length(matched_gene_ids),
                note = group_gene_filter_note(rnafilt_counts, gene_filter)
            )
        }
    }

    final_ids <- matching_group_sample_ids(
        clinic_annot = clinic_annot,
        rnafilt_counts = rnafilt_counts,
        clinic_filters = normalized_clinic_filters,
        group_gene_filter = normalized_gene_filters,
        sample_id_col = sample_id_col
    )
    final_ids <- intersect(available_ids, final_ids)

    list(
        available_ids = available_ids,
        clinic_filters = normalized_clinic_filters,
        group_gene_filter = normalized_gene_filters,
        clinic_ids = clinic_ids,
        gene_ids = gene_ids,
        final_ids = final_ids,
        gene_filter_counts = gene_filter_counts
    )
}


build_empty_group_error_message <- function(group_label, diagnostics) {
    lines <- c(
        sprintf("The %s-group filters do not match any sample after filtering.", group_label),
        sprintf("Samples available after global filtering: %d.", length(diagnostics$available_ids))
    )

    if (!is.null(diagnostics$clinic_filters)) {
        lines <- c(
            lines,
            sprintf(
                "Clinical rules matched %d sample(s): %s",
                length(diagnostics$clinic_ids),
                format_clinic_filter_description(diagnostics$clinic_filters)
            )
        )
    }

    if (!is.null(diagnostics$group_gene_filter)) {
        lines <- c(
            lines,
            sprintf(
                "Gene-expression rules matched %d sample(s) in combination: %s",
                length(diagnostics$gene_ids),
                format_group_gene_filter_description(diagnostics$group_gene_filter)
            )
        )

        per_gene_lines <- vapply(diagnostics$gene_filter_counts, function(entry) {
            line <- sprintf(
                "%s -> %d sample(s)",
                format_group_gene_filter_description(entry$filter),
                entry$match_count
            )

            if (!is.null(entry$note) && nzchar(entry$note)) {
                line <- sprintf("%s (%s)", line, entry$note)
            }

            line
        }, character(1))

        lines <- c(lines, sprintf("Per-gene rule matches: %s", paste(per_gene_lines, collapse = "; ")))
    }

    if (!is.null(diagnostics$clinic_filters) &&
        !is.null(diagnostics$group_gene_filter) &&
        length(diagnostics$clinic_ids) > 0 &&
        length(diagnostics$gene_ids) > 0 &&
        length(diagnostics$final_ids) == 0) {
        lines <- c(lines, "The clinical rules and gene-expression rules each match samples on their own, but their intersection is empty.")
    }

    paste(lines, collapse = "\n")
}


DEGSEA_pipe <- function(rnafilt_counts, clinic_annot, control_filters, test_filters, output_dir, pathways_to_use,
                        filter_by_gene = NULL, keep_low_or_high = NULL, quantile_thr = NULL, covariates = NULL,
                        clinic_filters = NULL, control_gene_filter = NULL, test_gene_filter = NULL, min_size = 15)
    {
    library(DESeq2)
    library(dplyr)
    library(GSVA)
    library(ggplot2)
    library(tidyverse)
    library(jsonlite)
    library(fgsea)
    incProgress(0.2)

    condition_col = "DEGSEA_group"
    sample_ids = colnames(rnafilt_counts)

    if (is.null(sample_ids)) {
        stop("rnafilt_counts must provide sample IDs in its column names.")
    }

    colnames(rnafilt_counts) = trimws(as.character(sample_ids))

    if (any(is.na(colnames(rnafilt_counts)) | !nzchar(colnames(rnafilt_counts)))) {
        stop("rnafilt_counts contains missing or empty sample IDs in its column names.")
    }

    control_filters = normalize_clinic_filters(control_filters)
    test_filters = normalize_clinic_filters(test_filters)
    control_gene_filter = normalize_group_gene_filters(control_gene_filter)
    test_gene_filter = normalize_group_gene_filters(test_gene_filter)

    if (is.null(control_filters) && is.null(control_gene_filter)) {
        stop("At least one control-group clinic or gene filter is required.")
    }

    if (is.null(test_filters) && is.null(test_gene_filter)) {
        stop("At least one test-group clinic or gene filter is required.")
    }

    comparison_suffix = comparison_filters_path_suffix(
        control_filters = control_filters,
        test_filters = test_filters,
        control_gene_filter = control_gene_filter,
        test_gene_filter = test_gene_filter
    )

    if(!is.null(covariates))
    {
        vars <- c(covariates, condition_col)
        DESeq_design <- reformulate(vars)
    }else{DESeq_design = reformulate(condition_col)
          vars = condition_col}


    pathways_name = pathways_to_use
    selected_pathways = normalize_gene_set_collection(get(pathways_to_use), fallback_prefix = pathways_to_use)

    if (length(selected_pathways) == 0) {
        stop(sprintf("The selected pathway collection '%s' does not contain any usable gene set.", pathways_to_use))
    }

    parsed_design = paste0(gsub("~", "", gsub(" ", "", DESeq_design))[2], "__", comparison_suffix)  # Will be used for the path of the save directory

    # Filter on clinic and gene selections before downstream analyses.
    filtered_data = filtering_on_clinic_and_genes(clinic_annot, rnafilt_counts, 
                                                  clinic_filters = clinic_filters,
                                                  filter_by_gene = filter_by_gene, 
                                                  keep_low_or_high = keep_low_or_high, 
                                                  quantile_thr = quantile_thr,
                                                  parsed_design = parsed_design,
                                                  output_dir =  output_dir,
                                                  folder_name = "DESeq") 
    
    rnafilt_counts = filtered_data$rnafilt_counts
    clinic_annot = ensure_clinic_sample_id_col(filtered_data$clinic_annot)
    output_DESeq = filtered_data$output_path
    filter_suffix = basename(dirname(output_DESeq))
    available_sample_ids = intersect(clinic_annot$ID_Patient, colnames(rnafilt_counts))

    if (length(available_sample_ids) == 0) {
        stop(
            "The global clinic/gene filtering step removed all samples before the control/test groups were defined.\n",
            "Adjust the filters in 'Filtering on clinic' or 'Filtering on gene'."
        )
    }
    
    incProgress(0.25)

    #### Remove NA and align clinic annot and RNAseq
    clinic_sample_ids = clinic_annot$ID_Patient
    control_ids = matching_group_sample_ids(
        clinic_annot = clinic_annot,
        rnafilt_counts = rnafilt_counts,
        clinic_filters = control_filters,
        group_gene_filter = control_gene_filter
    )
    test_ids = matching_group_sample_ids(
        clinic_annot = clinic_annot,
        rnafilt_counts = rnafilt_counts,
        clinic_filters = test_filters,
        group_gene_filter = test_gene_filter
    )

    control_ids = intersect(control_ids, colnames(rnafilt_counts))
    test_ids = intersect(test_ids, colnames(rnafilt_counts))

    if (length(control_ids) == 0) {
        stop(build_empty_group_error_message(
            "control",
            diagnose_group_matching(
                clinic_annot = clinic_annot,
                rnafilt_counts = rnafilt_counts,
                clinic_filters = control_filters,
                group_gene_filter = control_gene_filter
            )
        ))
    }

    if (length(test_ids) == 0) {
        stop(build_empty_group_error_message(
            "test",
            diagnose_group_matching(
                clinic_annot = clinic_annot,
                rnafilt_counts = rnafilt_counts,
                clinic_filters = test_filters,
                group_gene_filter = test_gene_filter
            )
        ))
    }

    overlapping_ids = intersect(control_ids, test_ids)
    if (length(overlapping_ids) > 0) {
        stop(sprintf(
            "Some samples match both control and test definitions: %s",
            paste(head(overlapping_ids, 10), collapse = ", ")
        ))
    }

    selected_ids = union(control_ids, test_ids)
    clinicannot_noNA = clinic_annot[
        clinic_sample_ids %in% selected_ids &
        clinic_sample_ids %in% colnames(rnafilt_counts),
        ,
        drop = FALSE
    ]
    clinicannot_sample_ids = clinicannot_noNA$ID_Patient
    clinicannot_noNA[[condition_col]] = ifelse(clinicannot_sample_ids %in% control_ids, "control", "test")

    if (!is.null(covariates)) {
        for (covar in covariates) {
            clinicannot_noNA = clinicannot_noNA[
                !is.na(clinicannot_noNA[[covar]]) &
                clinicannot_noNA$ID_Patient %in% colnames(rnafilt_counts),
            ]
        }
    }

    rnafilt_noNA <- rnafilt_counts[
    , colnames(rnafilt_counts) %in% clinicannot_noNA$ID_Patient,
    drop = FALSE
    ]

    rnafilt_noNA <- rnafilt_noNA[
    , rownames(clinicannot_noNA),
    drop = FALSE
    ]
    clinicannot_noNA[[condition_col]] = factor(clinicannot_noNA[[condition_col]], levels = c("control", "test"))

    control_ids = intersect(control_ids, rownames(clinicannot_noNA))
    test_ids = intersect(test_ids, rownames(clinicannot_noNA))

    if (length(control_ids) == 0 || length(test_ids) == 0) {
        stop("At least one DESeq group became empty after applying covariate filtering.")
    }

    incProgress(0.3)

    ## RUN DESeq ##
    DESeq_dds = doDGEv2(rnamat = round(rnafilt_noNA),
                      annot = clinicannot_noNA,
                      design = DESeq_design,
                      condition = condition_col,
                      modalities_control = "control",
                      modalities_test = "test"
                      )
                      
    # doDGEv2 collapses the selected modalities into the binary levels
    # "control" and "test" inside the DESeq2 object.
    DESeq_res = results(DESeq_dds, contrast = c(condition_col, "test", "control"))


    ## SAVE_RESULTS ##

    if (!dir.exists(output_DESeq)) dir.create(output_DESeq, recursive = TRUE)
    saveRDS(DESeq_dds, file = paste0(output_DESeq,"/DESeq_dds", ".rds"), compress = "xz")  # good compression

    incProgress(0.7)

    ## VOLCANO PLOT ## 
    log2FC_threshold = 0
    pval_threshold = 0.05

    df = data.frame(
    Gene = rownames(DESeq_res),
    log2FoldChange = DESeq_res$log2FoldChange,
    padj = DESeq_res$padj
    )

    tt_group = gsub("~", "", DESeq_design)[2]

    df$log2FoldChange[is.na(df$log2FoldChange)] = 0
    df$padj[is.na(df$padj)] = 1 

    df$color = "NS"  # Non significatif
    df$color[df$log2FoldChange <= -log2FC_threshold & df$padj < pval_threshold] = "Down"
    df$color[df$log2FoldChange >= log2FC_threshold & df$padj < pval_threshold] = "Up"

    df_down = df[which(df$log2FoldChange <0),]
    top_genes_down = df_down[order(df_down$padj, decreasing = FALSE), ][1:500, ] ### top 50 meilleurs gene auront leur label d'écrit
    df_up = df[which(df$log2FoldChange >0),]
    top_genes_up = df_up[order(df_up$padj, decreasing = FALSE), ][1:500, ] ### top 50 meilleurs gene auront leur label d'écrit
    top_genes= rbind(top_genes_down[1:50, ], top_genes_up[1:50, ])

    p = ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = color)) +
        geom_point(alpha = 0.4) +
        scale_color_manual(values = c("Down" = "#4ab3d6", "Up" = "lightcoral", "NS" = "grey")) +
        geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), linetype = "dashed") +
        geom_text(data = top_genes, aes(label = Gene), vjust = -0.5, size = 3) +
        theme_minimal() +
        labs(title = paste("Differential gene expression on", gsub("~", "", DESeq_design) ,"factor")[2], x = "Log2 Fold Change", y = "padj") +
        theme(legend.position = "none") +
        annotate(geom = 'text', label = "Up in control", x = -Inf, y = 0, hjust = 0, vjust = 0) +
        annotate(geom = 'text', label = "Up in test", x = Inf, y = 0, hjust = 1, vjust = 0)


    paste0(output_DESeq, "top_genes_down_500.csv")
    write_csv_mkdir(top_genes_down, paste0(output_DESeq, "/top_genes_down_500.csv"), row.names=FALSE)
    write_csv_mkdir(top_genes_up, paste0(output_DESeq, "/top_genes_up_500.csv"), row.names=FALSE)
    ggsave(paste0(output_DESeq, "/", "volcano_plot.jpg"), plot = p, width = 8, height = 8, dpi = 300)


    incProgress(0.75)

    #### GSEA ####

    output_GSEA = file.path(output_dir, "GSEA", filter_suffix, parsed_design)


    gene_stat = DESeq_res$stat
    names(gene_stat) = rownames(DESeq_res)
    ranked_gene_vec = sort(gene_stat, decreasing = TRUE)

    ranked_gene_df = as.data.frame(ranked_gene_vec, col.names = FALSE)
    ranked_gene_df$Genes = rownames(ranked_gene_df)
    colnames(ranked_gene_df) = c("Score", "Genes")
    print("GSEA classes are: control VS test")

    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    write.table(ranked_gene_df, paste0(output_DESeq, "/ranked_genes.rnk"))

    gsea_results = gsea_multi(ranked_gene_vec, selected_pathways, min_size = min_size)

    pvals = gsea_results$pval
    names(pvals) = gsea_results$pathway
    sortedpvals = sort(pvals)

    sorted_gsea_results = gsea_results[order(abs(gsea_results$NES), decreasing=TRUE), ]

    top_pathways = sorted_gsea_results[1:500, ]
    leading_genes = top_pathways$leadingEdge
    leading_genes = unique(unlist(strsplit(unlist(leading_genes), split = ",")))
    top_lead_genes = as.list(top_pathways[["leadingEdge"]])
    names(top_lead_genes) = top_pathways[["pathway"]]

    top_lead_genes = lapply(top_lead_genes, function(x) {
    unlist(strsplit(x, ",\\s*"))
    })
    write_csv_mkdir(sorted_gsea_results[, c("pathway", "NES", "pval", "padj", "log2err", "ES", "size", "leadingEdge")], paste0(output_GSEA, "/es_", pathways_name, ".csv"), row.names=TRUE)

    incProgress(0.95)

    #### GESECA ####

    output_GESECA = file.path(output_dir, "GESECA", filter_suffix, parsed_design, paste0(pathways_name, ".csv"))

    
    # On les met en no_NA pour la colonne de réponse
    rnafilt_norm = normVST_bulk(round(rnafilt_counts))
    clinicannot_noNA = clinic_annot[
        clinic_sample_ids %in% selected_ids &
        clinic_sample_ids %in% colnames(rnafilt_norm),
        ,
        drop = FALSE
    ]
    rnafilt_norm_noNA = rnafilt_norm[, colnames(rnafilt_norm) %in% clinicannot_noNA$ID_Patient]
    rnafilt_norm_noNA = rnafilt_norm_noNA[, rownames(clinicannot_noNA)]  # On les met dans le même ordre
    rnafilt_norm_cent_noNA = scale(rnafilt_norm_noNA)


    geseca_results = as.data.frame(geseca(selected_pathways, rnafilt_norm_cent_noNA, minSize = min_size, maxSize = 100000, eps = 0))
    rownames(geseca_results) = geseca_results$pathway

    ordered_ss = geseca_results[order(geseca_results$pctVar, decreasing = TRUE), ]

    top_geseca_names = ordered_ss$pathway[1:50]

    write_csv_mkdir(ordered_ss, output_GESECA, row.names=FALSE)



    #### ssGSEA ####

    output_ssGSEA = file.path(output_dir, "ssGSEA", filter_suffix, parsed_design, paste0("es_ssgsea_", pathways_name, ".csv"))

    sp = ssgseaParam(
    as.matrix(rnafilt_noNA),
    selected_pathways,
    minSize = min_size,
    alpha = 0.25,
    normalize = TRUE   # = active la normalisation ssGSEA
    )

    nes_ssgsea = gsva(sp)

    write_csv_mkdir(nes_ssgsea, output_ssGSEA, row.names = TRUE)

    return(list(
                plots = list(volcano = p),
                tables = list(top_DE = top_genes, 
                              top_GSEA = sorted_gsea_results[, c("pathway", "NES", "padj", "leadingEdge")], 
                              top_GESECA = ordered_ss, 
                              top_ssGSEA = nes_ssgsea)
                )
            )
    }
