set.seed(69)


adjust_on_ESTIMATE =  "immune_score+stromal_score"  # = stromal_score, immune_score, estimate_score


anchored_GSEA_pipe <- function(rnafilt_counts, clinic_annot, output_dir, anchor_gene, pathways_to_use,
                        filter_by_gene, keep_low_or_high, quantile_thr = NULL, adjust_on_ESTIMATE = NULL,
                        clinic_filters = NULL, min_size = 10)
    {
    library(data.table)
    library(fgsea)
    library(tidyestimate)
    library(limma)

    adjustment_suffix = gsub("_score", "", adjust_on_ESTIMATE)
    selected_pathways = normalize_gene_set_collection(get(pathways_to_use), fallback_prefix = pathways_to_use)

    if (length(selected_pathways) == 0) {
        stop(sprintf("The selected pathway collection '%s' does not contain any usable gene set.", pathways_to_use))
    }
    # Load datasets

    # Filter clinic annotations
    # clinic_annot = clinic_annot[clinic_annot[, subset_column] %in% subset_modality,]
    # rnafilt_counts = rnafilt_counts[, colnames(rnafilt_counts) %in% rownames(clinic_annot)]
    filtered_data = filtering_on_clinic_and_genes(clinic_annot, rnafilt_counts, 
                                  clinic_filters = clinic_filters,
                                  filter_by_gene = filter_by_gene, 
                                  keep_low_or_high = keep_low_or_high, 
                                  quantile_thr = quantile_thr, 
                                  parsed_design = adjustment_suffix,
                                  output_dir = output_dir,
                                  folder_name = "anchored_GSEA") 
    
    rnafilt_counts = filtered_data$rnafilt_counts
    clinic_annot = filtered_data$clinic_annot
    output_path = paste0(filtered_data$output_path, "_", anchor_gene)

    if (!is.null(adjust_on_ESTIMATE) && length(adjust_on_ESTIMATE) > 0) {
    # TPM: genes x samples, rownames = HGNC symbols, values = TPM (NOT log)
        df <- data.frame(gene = rownames(rnafilt_counts), rnafilt_counts, check.names = FALSE)
        est <- estimate_score(df, is_affymetrix = FALSE)  # immune/stromal/estimate scores

        E <- rnafilt_counts
        meta <- clinic_annot[colnames(E), ]       
        meta[[anchor_gene]] <- as.numeric(E[anchor_gene, ])
        meta$immune_score = as.numeric(est$immune)
        meta$stromal_score = as.numeric(est$stromal)
        meta$estimate_score = as.numeric(est$estimate)

        # Put whatever you have here: cohort/batch + purity + cell fractions
        formula <- as.formula(paste0("~ scale(", anchor_gene, ") + ", adjust_on_ESTIMATE))
        design <- model.matrix(formula, data = meta)

        fit <- eBayes(lmFit(E, design), robust=TRUE)
        tt <- topTable(fit, coef=paste0("scale(",anchor_gene,")"), number=Inf, sort.by="t")
        ranks <- tt$t
        names(ranks) <- rownames(tt)
        ranks <- ranks[names(ranks) != anchor_gene]  # On enleve BSG prck forcement il est corrélé à lui même

        res_adjusted <- fgseaMultilevel(pathways=selected_pathways, stats=ranks, minSize=min_size, maxSize=10000)
        res_adjusted <- res_adjusted[order(abs(res_adjusted$padj), decreasing = FALSE), ]
        res_adjusted = flatten_list_cols(res_adjusted)

        write_csv_mkdir(res_adjusted, paste0(output_path, ".csv"))
        
        return(res_adjusted)

    }else {
        # Genes are ranked by expression correlation with BSG 
        # NES < 0 : anticorrelated with BSG
        # NES > 0 : correlated with BSG
        E <- rnafilt_counts   # genes x samples, rownames = symbols
        anch_gene <- as.numeric(E[anchor_gene, ])

        # gene-wise association to BSG
        stats <- apply(E, 1, function(x)
        cor(x, anch_gene, method = "spearman", use = "pairwise.complete.obs")
        )

        stats <- sort(stats, decreasing = TRUE)
        stats <- stats[names(stats) != anchor_gene]  # On enleve BSG prck forcement il est corrélé à lui même

        # selected_pathways: named list of character vectors (gene symbols)
        res_simple <- fgseaMultilevel(pathways = selected_pathways, stats = stats, minSize = min_size, maxSize = 500)
        res_simple <- res_simple[order(abs(res_simple$padj), decreasing = FALSE), ]

        res_simple = flatten_list_cols(res_simple)

        write_csv_mkdir(res_simple, paste0(output_path, ".csv"))
        return(res_simple)

    }



    # out <- copy(fg_bsg_adj)
    # out <- out[order(-abs(NES))]

    # list_cols <- names(out)[vapply(out, is.list, logical(1))]

    # out[, (list_cols) := lapply(.SD, function(x)
    #   vapply(x, function(xx) paste(xx, collapse = ";"), character(1))
    # ), .SDcols = list_cols]

    # write_csv_mkdir()

    # fwrite(out,
    #        "/home/quentin/01_PROJETS/01_Predimel/07_BSG/BSG_GSEA/Filt_ES_42.csv",
    #        sep = ",", dec = ".")

}


