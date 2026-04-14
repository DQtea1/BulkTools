library(DESeq2)
library(dplyr)
library(GSVA)
library(ggplot2)
library(tidyverse)
library(jsonlite)
library(fgsea)


DEGSEA_pipe <- function(rnafilt_counts, clinic_annot, contrast, control, test, output_dir, pathways_to_use,
                        filter_by_gene = NULL, keep_low_or_high = NULL, quantile_thr = NULL, covariates = NULL, 
                        clinic_filters = NULL)
    {
    if(!is.null(covariates))
    {
        str(covariates)
        vars <- c(covariates, contrast)
        DESeq_design <- reformulate(vars)
    }else{DESeq_design = reformulate(contrast)
          vars = contrast}


    pathways_name = pathways_to_use
    selected_pathways = get(pathways_to_use) # "all_pathways" "REACTOME_pathways", "GOBP_pathways", "KEGG_pathways", "nerve_sigs", "nerve_signatures", "final_nerve_signatures"

    parsed_design = gsub("~", "", gsub(" ", "", DESeq_design))[2]  # Will be used for the path of the save directory

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
    clinic_annot = filtered_data$clinic_annot
    output_DESeq = filtered_data$output_path
    filter_suffix = basename(dirname(output_DESeq))

    #### Remove NA and align clinic annot and RNAseq
    rownames(clinic_annot) = clinic_annot$ID_Patient
    clinicannot_noNA = clinic_annot[!is.na(clinic_annot[[contrast]]) &
                                    clinic_annot$ID_Patient %in% colnames(rnafilt_counts), ]

    clinicannot_noNA = clinicannot_noNA[clinicannot_noNA[[contrast]] %in% c(control, test), ]

    if (!is.null(covariates)) {
        for (covar in vars) {
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
    clinicannot_noNA[[contrast]] = factor(clinicannot_noNA[[contrast]], levels = c(control, test))  


    ## RUN DESeq ##
    DESeq_dds = doDGEv2(rnamat = round(rnafilt_noNA),
                      annot = clinicannot_noNA,
                      design = DESeq_design,
                      condition = contrast,
                      modalities_control = control,
                      modalities_test = test
                      )
                      
    # doDGEv2 collapses the selected modalities into the binary levels
    # "control" and "test" inside the DESeq2 object.
    DESeq_res = results(DESeq_dds, contrast = c(contrast, "test", "control"))


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

    ctrl_group = "control"
    tt_group = gsub("~", "", DESeq_design)[2]

    df$log2FoldChange[is.na(df$log2FoldChange)] = 0
    df$padj[is.na(df$padj)] = 1 

    df$color = "NS"  # Non significatif
    df$color[df$log2FoldChange <= -log2FC_threshold & df$padj < pval_threshold] = "Down"
    df$color[df$log2FoldChange >= log2FC_threshold & df$padj < pval_threshold] = "Up"

    # df$direction = ifelse(df$log2FoldChange > 1, control, 
    #                            ifelse(df$log2FoldChange < -1, test, "NS"))

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
        annotate(geom = 'text', label = paste0("Up in ", control), x = -Inf, y = 0, hjust = 0, vjust = 0) +
        annotate(geom = 'text', label = paste0("Up in ", test), x = Inf, y = 0, hjust = 1, vjust = 0)


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
    print(paste("GSEA classes are:", control, "VS", test))

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
    clinicannot_noNA = clinic_annot[!is.na(clinic_annot[[contrast]])
                                        & clinic_annot$ID_Patient %in% colnames(rnafilt_norm), ]
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
