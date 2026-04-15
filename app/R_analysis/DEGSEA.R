library(DESeq2)
library(dplyr)
library(GSVA)
library(ggplot2)
library(tidyverse)
library(jsonlite)
library(fgsea)


DEGSEA_pipe <- function(rnafilt_counts, clinic_annot, control_filters, test_filters, output_dir, pathways_to_use,
                        filter_by_gene = NULL, keep_low_or_high = NULL, quantile_thr = NULL, covariates = NULL, 
                        clinic_filters = NULL, control_gene_filter = NULL, test_gene_filter = NULL)
    {
    condition_col = "DEGSEA_group"
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
        stop("The control-group filters do not match any sample after filtering.")
    }

    if (length(test_ids) == 0) {
        stop("The test-group filters do not match any sample after filtering.")
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

    gsea_results = gsea_multi(ranked_gene_vec, selected_pathways)

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


    geseca_results = as.data.frame(geseca(selected_pathways, rnafilt_norm_cent_noNA, minSize = 15, maxSize = 100000, eps = 0))
    rownames(geseca_results) = geseca_results$pathway

    ordered_ss = geseca_results[order(geseca_results$pctVar, decreasing = TRUE), ]

    top_geseca_names = ordered_ss$pathway[1:50]

    write_csv_mkdir(ordered_ss, output_GESECA, row.names=FALSE)



    #### ssGSEA ####

    output_ssGSEA = file.path(output_dir, "ssGSEA", filter_suffix, parsed_design, paste0("es_ssgsea_", pathways_name, ".csv"))

    sp = ssgseaParam(
    as.matrix(rnafilt_noNA),
    selected_pathways,
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
