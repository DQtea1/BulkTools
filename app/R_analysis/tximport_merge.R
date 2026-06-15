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
                                remove_gene_classes = NULL,
                                progress_cb = NULL) {
  library(tximport)
  library(ggplot2)

  # Optional progress reporter: progress_cb(frac, detail) with frac in [0, 1].
  # No-op when called outside Shiny (progress_cb left NULL).
  report <- function(frac, detail) {
    if (is.function(progress_cb)) progress_cb(frac, detail)
  }

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
    # Several suffixes can be given, separated by "|" (e.g. ".tsv|.txt|.csv").
    strip_tokens <- trimws(strsplit(file_strip, "|", fixed = TRUE)[[1]])
    strip_tokens <- strip_tokens[nzchar(strip_tokens)]
    # Strip the longest/most specific tokens first so e.g. "_quantif.tsv"
    # wins over ".tsv" when both are provided.
    strip_tokens <- strip_tokens[order(nchar(strip_tokens), decreasing = TRUE)]

    if (length(strip_tokens) == 0) {
      stop("'File suffix to strip' is non-empty but contains no usable token after splitting on '|'.")
    }

    keep <- Reduce(`|`, lapply(strip_tokens, function(tok) grepl(tok, bulk_filenames, fixed = TRUE)))
    if (!any(keep)) {
      stop(sprintf("No files matching any suffix [%s] in '%s'.",
                   paste(strip_tokens, collapse = ", "), bulk_folder))
    }
    bulk_filenames <- bulk_filenames[keep]

    sample_names <- bulk_filenames
    for (tok in strip_tokens) {
      sample_names <- sub(tok, "", sample_names, fixed = TRUE)
    }
  } else {
    sample_names <- tools::file_path_sans_ext(bulk_filenames)
  }
  bulk_filepaths <- file.path(bulk_folder, bulk_filenames)

  samples <- data.frame(ID_Patient = sample_names,
                        filepath   = bulk_filepaths,
                        stringsAsFactors = FALSE)
  bulk_files <- samples$filepath
  names(bulk_files) <- samples$ID_Patient

  ## 2) tximport per file, then merge at the gene level.
  ##    tximport's batch mode requires identical tx IDs across files; running
  ##    per file is robust to heterogeneous quant outputs (different reference,
  ##    missing transcripts) at the cost of one tximport call per sample.
  counts_list <- vector("list", length(bulk_files))
  tpm_list    <- vector("list", length(bulk_files))
  names(counts_list) <- names(bulk_files)
  names(tpm_list)    <- names(bulk_files)

  for (idx in seq_along(bulk_files)) {
    sid <- names(bulk_files)[idx]
    fp  <- bulk_files[[idx]]
    report(0.1 + 0.6 * (idx / length(bulk_files)),
           sprintf("tximport %d/%d : %s", idx, length(bulk_files), sid))
    ex  <- read.delim(fp, sep = "\t")
    if (ncol(ex) < 2) {
      stop(sprintf("Quantification file '%s' has fewer than 2 columns; cannot build tx2gene.", fp))
    }
    res_i <- tximport(
      setNames(fp, sid),
      type                = tx_type,
      txIn                = tx_txIn,
      countsFromAbundance = tx_countsFromAbundance,
      tx2gene             = ex[, 1:2],
      geneIdCol           = tx_geneIdCol,
      txIdCol             = tx_txIdCol,
      abundanceCol        = tx_abundanceCol,
      countsCol           = tx_countsCol,
      lengthCol           = tx_lengthCol,
      ignoreTxVersion     = tx_ignoreTxVersion
    )
    counts_list[[sid]] <- setNames(as.numeric(res_i$counts[, 1]),    rownames(res_i$counts))
    tpm_list[[sid]]    <- setNames(as.numeric(res_i$abundance[, 1]), rownames(res_i$abundance))
  }

  # Drop transcripts that summarized to a missing/empty gene name (NA or ""),
  # and collapse any duplicate gene names within a sample by summing. Without
  # this, an NA gene name survives in a sample vector but is dropped from the
  # sorted gene universe, causing a 'subscript out of bounds' on assignment.
  clean_named_vec <- function(x) {
    nm   <- names(x)
    keep <- !is.na(nm) & nzchar(nm)
    x    <- x[keep]
    if (anyDuplicated(names(x))) {
      agg <- tapply(x, names(x), sum)
      x   <- setNames(as.numeric(agg), names(agg))
    }
    x
  }
  counts_list <- lapply(counts_list, clean_named_vec)
  tpm_list    <- lapply(tpm_list,    clean_named_vec)

  all_genes <- sort(unique(Reduce(union, lapply(counts_list, names))))

  assemble_matrix <- function(vec_list, gene_universe, fill = 0) {
    m <- matrix(fill,
                nrow = length(gene_universe),
                ncol = length(vec_list),
                dimnames = list(gene_universe, names(vec_list)))
    for (sid in names(vec_list)) {
      v    <- vec_list[[sid]]
      nm   <- names(v)
      keep <- !is.na(nm) & nm %in% gene_universe
      if (any(keep)) m[nm[keep], sid] <- v[keep]
    }
    m
  }

  RNAseq_counts <- assemble_matrix(counts_list, all_genes, fill = 0)
  RNAseq_TPM    <- assemble_matrix(tpm_list,    all_genes, fill = 0)

  ## 3) Save the unfiltered triple
  unfiltered_dir <- file.path(output_dir, "00_Unfiltered")
  if (!dir.exists(unfiltered_dir)) dir.create(unfiltered_dir, recursive = TRUE)

  report(0.74, "VST on unfiltered matrix + writing")
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
  report(0.9, "VST on filtered matrix + writing")
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
