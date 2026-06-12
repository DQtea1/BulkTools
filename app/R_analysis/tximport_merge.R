# Wraps the bulk RNA-seq merge step:
# tximport per-sample quant files -> raw counts + TPM
# -> drop low-depth samples, drop low-count genes, optionally drop gene classes
# -> write unfiltered and filtered triples (counts / TPM / VST) to disk.

tximport_merge_pipe <- function(bulk_folder,
                                output_dir,
                                file_strip = "_quantif.tsv",
                                tx_type = "salmon",
                                tx_txIn = TRUE,
                                tx_countsFromAbundance = "lengthScaledTPM",
                                tx_geneIdCol = "Gene_name",
                                tx_txIdCol = "Name",
                                tx_abundanceCol = "TPM",
                                tx_countsCol = "NumReads",
                                tx_lengthCol = "Length",
                                tx_ignoreTxVersion = FALSE,
                                min_gene_count = 15,
                                min_n_samples = 0.33,
                                min_seq_depth = 30000000,
                                remove_gene_classes = NULL) {
  library(tximport)
  library(ggplot2)

  if (is.null(bulk_folder) || !nzchar(bulk_folder) || !dir.exists(bulk_folder)) {
    stop(sprintf("Bulk folder does not exist: '%s'", bulk_folder))
  }
  if (is.null(output_dir) || !nzchar(output_dir)) {
    stop("Output directory must be provided.")
  }

  ## 1) Discover input quant files
  bulk_filenames <- list.files(path = bulk_folder, full.names = FALSE)
  if (length(bulk_filenames) == 0) {
    stop(sprintf("No files found in '%s'.", bulk_folder))
  }

  if (!is.null(file_strip) && nzchar(file_strip)) {
    keep <- grepl(file_strip, bulk_filenames, fixed = TRUE)
    if (!any(keep)) {
      stop(sprintf("No files matching suffix '%s' in '%s'.", file_strip, bulk_folder))
    }
    bulk_filenames <- bulk_filenames[keep]
    sample_names   <- sub(file_strip, "", bulk_filenames, fixed = TRUE)
  } else {
    sample_names <- tools::file_path_sans_ext(bulk_filenames)
  }
  bulk_filepaths <- file.path(bulk_folder, bulk_filenames)

  samples <- data.frame(ID_Patient = sample_names,
                        filepath   = bulk_filepaths,
                        stringsAsFactors = FALSE)
  bulk_files <- samples$filepath
  names(bulk_files) <- samples$ID_Patient

  ## 2) tximport
  exemple_bulk <- read.delim(bulk_filepaths[1], sep = "\t")
  if (ncol(exemple_bulk) < 2) {
    stop("Quantification files have fewer than 2 columns; cannot derive tx2gene from the first two columns.")
  }

  RNAseq_merged <- tximport(
    bulk_files,
    type                = tx_type,
    txIn                = tx_txIn,
    countsFromAbundance = tx_countsFromAbundance,
    tx2gene             = exemple_bulk[, 1:2],
    geneIdCol           = tx_geneIdCol,
    txIdCol             = tx_txIdCol,
    abundanceCol        = tx_abundanceCol,
    countsCol           = tx_countsCol,
    lengthCol           = tx_lengthCol,
    ignoreTxVersion     = tx_ignoreTxVersion
  )

  RNAseq_counts <- RNAseq_merged$counts
  RNAseq_TPM    <- RNAseq_merged$abundance

  ## 3) Save the unfiltered triple
  unfiltered_dir <- file.path(output_dir, "00_Unfiltered")
  if (!dir.exists(unfiltered_dir)) dir.create(unfiltered_dir, recursive = TRUE)

  RNAseq_unfilt_vst <- normVST_bulk(round(RNAseq_counts))

  write_csv_mkdir(RNAseq_counts,     file.path(unfiltered_dir, "RNAseq_counts.csv"))
  write_csv_mkdir(RNAseq_TPM,        file.path(unfiltered_dir, "RNAseq_TPM.csv"))
  write_csv_mkdir(RNAseq_unfilt_vst, file.path(unfiltered_dir, "RNAseq_vst.csv"))

  ## 4) Per-sample read depth + plot
  total_counts <- colSums(as.matrix(RNAseq_counts), na.rm = TRUE)
  total_counts_df <- data.frame(
    ID_Patient       = names(total_counts),
    total_read_count = as.numeric(total_counts),
    stringsAsFactors = FALSE,
    row.names        = NULL
  )
  total_counts_df <- total_counts_df[order(total_counts_df$total_read_count), ]
  total_counts_df$kept <- total_counts_df$total_read_count >= min_seq_depth
  total_counts_df$ID_Patient_fac <- factor(total_counts_df$ID_Patient,
                                           levels = total_counts_df$ID_Patient)

  depth_plot <- ggplot(total_counts_df,
                       aes(x = ID_Patient_fac, y = total_read_count, fill = kept)) +
    geom_col() +
    geom_hline(yintercept = min_seq_depth, linetype = "dashed", color = "black") +
    scale_fill_manual(values = c(`TRUE` = "#4caf50", `FALSE` = "#e57373"),
                      labels = c(`TRUE` = "kept", `FALSE` = "dropped"),
                      name = NULL) +
    labs(title = "Total read depth per sample",
         subtitle = sprintf("dashed line = min_seq_depth = %s",
                            format(min_seq_depth, big.mark = ",", scientific = FALSE)),
         x = NULL, y = "Total read count") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))

  ggsave(file.path(unfiltered_dir, "read_depth_per_sample.png"),
         depth_plot, width = 14, height = 6, dpi = 200)

  kept_patients    <- total_counts_df$ID_Patient[ total_counts_df$kept]
  dropped_patients <- total_counts_df$ID_Patient[!total_counts_df$kept]

  if (length(kept_patients) < 2) {
    stop(sprintf(
      "Sequencing-depth filtering left %d sample(s); need at least 2 to continue (VST will fail).",
      length(kept_patients)
    ))
  }

  RNAseq_depth_counts <- RNAseq_counts[, kept_patients, drop = FALSE]
  RNAseq_depth_TPM    <- RNAseq_TPM[,    kept_patients, drop = FALSE]

  ## 5) Low-count gene filter
  M <- as.matrix(RNAseq_depth_counts)
  N <- ncol(M) * min_n_samples
  keep_genes <- rowSums(M >= min_gene_count) >= N
  RNAseq_filt_counts <- M[keep_genes, , drop = FALSE]
  n_genes_after_count_filter <- nrow(RNAseq_filt_counts)

  ## 6) Drop gene classes (RPS / RPL / MT) before normalization
  removed_classes <- character(0)
  if (!is.null(remove_gene_classes) && length(remove_gene_classes) > 0) {
    class_patterns <- list(RPS = "^RPS", RPL = "^RPL", MT = "^MT-")
    for (cls in remove_gene_classes) {
      pat <- class_patterns[[cls]]
      if (is.null(pat)) next
      n_before <- nrow(RNAseq_filt_counts)
      RNAseq_filt_counts <- RNAseq_filt_counts[!grepl(pat, rownames(RNAseq_filt_counts)), , drop = FALSE]
      removed_classes <- c(removed_classes,
                           sprintf("%s (-%d)", cls, n_before - nrow(RNAseq_filt_counts)))
    }
  }

  RNAseq_filt_TPM <- RNAseq_depth_TPM[rownames(RNAseq_filt_counts), , drop = FALSE]
  RNAseq_filt_vst <- normVST_bulk(round(RNAseq_filt_counts))

  ## 7) Save the filtered triple
  filtered_dir <- file.path(
    output_dir,
    sprintf("01_Filtered_min%s_in%s", min_gene_count, min_n_samples)
  )
  if (!dir.exists(filtered_dir)) dir.create(filtered_dir, recursive = TRUE)

  write_csv_mkdir(RNAseq_filt_counts, file.path(filtered_dir, "RNAseq_counts.csv"))
  write_csv_mkdir(RNAseq_filt_TPM,    file.path(filtered_dir, "RNAseq_TPM.csv"))
  write_csv_mkdir(RNAseq_filt_vst,    file.path(filtered_dir, "RNAseq_vst.csv"))

  ## 8) Build summary text
  dropped_preview <- if (length(dropped_patients) == 0) {
    "none"
  } else if (length(dropped_patients) <= 20) {
    paste(dropped_patients, collapse = ", ")
  } else {
    sprintf("%s, ... (+%d more)",
            paste(head(dropped_patients, 20), collapse = ", "),
            length(dropped_patients) - 20)
  }

  summary_text <- paste(
    sprintf("Input files matching '%s' : %d", file_strip, length(bulk_files)),
    sprintf("Unfiltered matrix          : %d genes x %d samples",
            nrow(RNAseq_counts), ncol(RNAseq_counts)),
    "",
    sprintf("Sequencing depth filter   (>= %s reads): kept %d / %d samples",
            format(min_seq_depth, big.mark = ",", scientific = FALSE),
            length(kept_patients), ncol(RNAseq_counts)),
    sprintf("  dropped samples         : %s", dropped_preview),
    "",
    sprintf("Gene count filter         (>= %d reads in >= %.0f%% of samples): kept %d genes",
            min_gene_count, min_n_samples * 100, n_genes_after_count_filter),
    if (length(removed_classes) > 0)
      sprintf("Gene class drops          : %s", paste(removed_classes, collapse = "; "))
    else
      "Gene class drops          : none",
    "",
    sprintf("Filtered matrix           : %d genes x %d samples",
            nrow(RNAseq_filt_counts), ncol(RNAseq_filt_counts)),
    "",
    "Files written:",
    sprintf("  %s/{RNAseq_counts,RNAseq_TPM,RNAseq_vst}.csv", unfiltered_dir),
    sprintf("  %s/{RNAseq_counts,RNAseq_TPM,RNAseq_vst}.csv", filtered_dir),
    sep = "\n"
  )

  return(list(
    plots   = list(depth = depth_plot),
    summary = summary_text,
    tables  = list(
      unfiltered_counts = RNAseq_counts,
      filtered_counts   = RNAseq_filt_counts,
      depth_per_sample  = total_counts_df
    )
  ))
}
