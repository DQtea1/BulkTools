# R/helpers_io.R
# Contains :
# read_delim_auto()
# write_csv_mkdir()
# read_tables_from_paths()
# filtering_on_clinic_and_genes()

read_delim_auto <- function(fileinfo) {
  shiny::req(fileinfo)
  ext <- tolower(tools::file_ext(fileinfo$name))
  if (ext == "csv") {
    read.csv(fileinfo$datapath, row.names = 1, check.names = FALSE)
  } else if (ext %in% c("tsv", "txt")) {
    read.delim(fileinfo$datapath, row.names = 1, check.names = FALSE)
  } else {
    stop("Unsupported file type: ", ext)
  }
}


write_csv_mkdir <- function(x, file, ...) {
  # If 'file' is a character path, ensure its parent dir exists
  if (is.character(file)) {
    dir <- dirname(file)
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  }
  write.csv(x, file = file, ...)
  invisible(file)
}


extract_clinic_sample_ids <- function(clinic_annot, sample_id_col = "ID_Patient") {
  if (sample_id_col %in% colnames(clinic_annot)) {
    sample_ids <- clinic_annot[[sample_id_col]]
  } else if (!is.null(rownames(clinic_annot))) {
    sample_ids <- rownames(clinic_annot)
  } else {
    stop(sprintf(
      "Sample IDs were not found in clinic_annot: column '%s' is missing and row names are empty.",
      sample_id_col
    ))
  }

  sample_ids <- trimws(as.character(sample_ids))

  if (length(sample_ids) != nrow(clinic_annot)) {
    stop("The extracted clinic sample IDs do not match the number of rows in clinic_annot.")
  }

  if (any(is.na(sample_ids) | !nzchar(sample_ids))) {
    stop("clinic_annot contains missing or empty sample IDs.")
  }

  sample_ids
}


ensure_clinic_sample_id_col <- function(clinic_annot, sample_id_col = "ID_Patient") {
  sample_ids <- extract_clinic_sample_ids(clinic_annot, sample_id_col = sample_id_col)
  clinic_annot[[sample_id_col]] <- sample_ids
  rownames(clinic_annot) <- sample_ids
  clinic_annot
}


normalize_clinic_filters <- function(clinic_filters) {
  if (is.null(clinic_filters) || length(clinic_filters) == 0) {
    return(NULL)
  }

  filter_names <- names(clinic_filters)
  if (is.null(filter_names)) {
    stop("clinic_filters must be a named list.")
  }

  normalized_filters <- list()

  for (idx in seq_along(clinic_filters)) {
    filter_name <- trimws(as.character(filter_names[[idx]]))

    if (!nzchar(filter_name)) {
      next
    }

    raw_values <- clinic_filters[[idx]]
    values <- as.character(unlist(raw_values, recursive = TRUE, use.names = FALSE))
    values <- trimws(values)
    values <- values[!is.na(values) & nzchar(values)]

    if (length(values) == 0) {
      next
    }

    existing_values <- normalized_filters[[filter_name]]
    normalized_filters[[filter_name]] <- unique(c(existing_values, values))
  }

  if (length(normalized_filters) == 0) {
    return(NULL)
  }

  normalized_filters[sort(names(normalized_filters))]
}


sanitize_filter_token <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("[/\\\\:*?\"<>|]+", "-", x)
  x <- gsub("\\s+", "_", x)
  x <- gsub("-{2,}", "-", x)
  x
}


stable_path_hash <- function(x) {
  raw_values <- utf8ToInt(enc2utf8(paste(x, collapse = "|")))

  if (length(raw_values) == 0) {
    return("000000000000")
  }

  mod_1 <- 2147483629
  mod_2 <- 2147483587
  hash_1 <- 0
  hash_2 <- 7

  for (idx in seq_along(raw_values)) {
    value <- raw_values[[idx]]
    hash_1 <- (hash_1 * 131 + value + idx) %% mod_1
    hash_2 <- (hash_2 * 137 + value + 3 * idx) %% mod_2
  }

  paste0(
    sprintf("%06x", as.integer(hash_1 %% 16777215)),
    sprintf("%06x", as.integer(hash_2 %% 16777215))
  )
}


compact_path_component <- function(x, max_chars = 140, prefix_chars = 96) {
  if (is.null(x) || length(x) == 0) {
    return("undefined")
  }

  x <- paste(x, collapse = "__")
  x <- sanitize_filter_token(x)

  if (nchar(x, type = "chars") <= max_chars) {
    return(x)
  }

  prefix_chars <- min(prefix_chars, max_chars - 16)
  paste0(
    substr(x, 1, prefix_chars),
    "__h",
    stable_path_hash(x)
  )
}


compact_filter_values <- function(values, preview_values = 3, max_chars = 80) {
  if (length(values) == 0) {
    return("empty")
  }

  values <- sanitize_filter_token(sort(unique(values)))
  joined_values <- paste(values, collapse = "+")

  if (length(values) <= preview_values && nchar(joined_values, type = "chars") <= max_chars) {
    return(joined_values)
  }

  shown_values <- values[seq_len(min(preview_values, length(values)))]

  paste0(
    paste(shown_values, collapse = "+"),
    "+more",
    length(values) - length(shown_values),
    "-h",
    substr(stable_path_hash(joined_values), 1, 8)
  )
}


clinic_filters_path_suffix <- function(clinic_filters) {
  normalized_filters <- normalize_clinic_filters(clinic_filters)

  if (is.null(normalized_filters)) {
    return("nofilt")
  }

  suffix_parts <- vapply(names(normalized_filters), function(filter_name) {
    values <- sort(unique(normalized_filters[[filter_name]]))
    compact_path_component(paste0(
      sanitize_filter_token(filter_name),
      "-",
      compact_filter_values(values)
    ), max_chars = 96, prefix_chars = 72)
  }, character(1))

  compact_path_component(paste(suffix_parts, collapse = "__"), max_chars = 120, prefix_chars = 88)
}


normalize_group_gene_filter <- function(group_gene_filter) {
  if (is.null(group_gene_filter)) {
    return(NULL)
  }

  gene_name <- trimws(as.character(group_gene_filter$gene))
  if (length(gene_name) == 0 || is.na(gene_name) || !nzchar(gene_name)) {
    return(NULL)
  }

  quantile_thr <- suppressWarnings(as.numeric(group_gene_filter$quantile_thr)[1])
  keep_low_or_high <- trimws(as.character(group_gene_filter$keep_low_or_high)[1])

  if (is.na(quantile_thr) || quantile_thr < 0 || quantile_thr > 1) {
    stop("Group gene-filter quantile must be a single numeric value between 0 and 1.")
  }

  if (!keep_low_or_high %in% c("low", "high")) {
    stop("Group gene-filter direction must be either 'low' or 'high'.")
  }

  list(
    gene = gene_name,
    quantile_thr = quantile_thr,
    keep_low_or_high = keep_low_or_high
  )
}


normalize_group_gene_filters <- function(group_gene_filters) {
  if (is.null(group_gene_filters) || length(group_gene_filters) == 0) {
    return(NULL)
  }

  is_single_filter <- is.list(group_gene_filters) &&
    !is.null(names(group_gene_filters)) &&
    any(names(group_gene_filters) %in% c("gene", "quantile_thr", "keep_low_or_high"))

  raw_filters <- if (is_single_filter) {
    list(group_gene_filters)
  } else {
    group_gene_filters
  }

  normalized_filters <- lapply(raw_filters, normalize_group_gene_filter)
  normalized_filters <- Filter(Negate(is.null), normalized_filters)

  if (length(normalized_filters) == 0) {
    return(NULL)
  }

  normalized_filters
}


gene_filter_path_suffix <- function(group_gene_filter) {
  normalized_gene_filter <- normalize_group_gene_filter(group_gene_filter)

  if (is.null(normalized_gene_filter)) {
    return(NULL)
  }

  compact_path_component(paste(
    "gene",
    sanitize_filter_token(normalized_gene_filter$gene),
    sanitize_filter_token(normalized_gene_filter$keep_low_or_high),
    paste0("q", sanitize_filter_token(normalized_gene_filter$quantile_thr)),
    sep = "-"
  ), max_chars = 72, prefix_chars = 56)
}


group_gene_filters_path_suffix <- function(group_gene_filters) {
  normalized_filters <- normalize_group_gene_filters(group_gene_filters)

  if (is.null(normalized_filters)) {
    return(NULL)
  }

  suffix_parts <- vapply(normalized_filters, gene_filter_path_suffix, character(1))
  suffix_parts <- suffix_parts[!is.na(suffix_parts) & nzchar(suffix_parts)]

  if (length(suffix_parts) == 0) {
    return(NULL)
  }

  compact_path_component(paste(suffix_parts, collapse = "__"), max_chars = 120, prefix_chars = 88)
}


matching_clinic_sample_ids <- function(clinic_annot, clinic_filters, sample_id_col = "ID_Patient") {
  normalized_filters <- normalize_clinic_filters(clinic_filters)

  if (is.null(normalized_filters)) {
    return(character(0))
  }
  sample_ids <- extract_clinic_sample_ids(clinic_annot, sample_id_col = sample_id_col)

  keep_samples <- rep(TRUE, nrow(clinic_annot))

  for (filter_name in names(normalized_filters)) {
    if (!filter_name %in% colnames(clinic_annot)) {
      stop(sprintf("Clinical column '%s' not found in clinic_annot.", filter_name))
    }

    keep_samples <- keep_samples & clinic_annot[, filter_name] %in% normalized_filters[[filter_name]]
  }

  unique(sample_ids[keep_samples])
}


matching_gene_sample_ids <- function(rnafilt_counts, group_gene_filter) {
  normalized_gene_filter <- normalize_group_gene_filter(group_gene_filter)

  if (is.null(normalized_gene_filter)) {
    return(character(0))
  }

  if (is.data.frame(rnafilt_counts)) {
    rnafilt_counts <- as.matrix(rnafilt_counts)
  }

  if (!normalized_gene_filter$gene %in% rownames(rnafilt_counts)) {
    stop(sprintf("Gene '%s' not found in rnafilt_counts.", normalized_gene_filter$gene))
  }

  gene_expr <- as.numeric(rnafilt_counts[normalized_gene_filter$gene, ])
  thr <- quantile(gene_expr, probs = normalized_gene_filter$quantile_thr, na.rm = TRUE)

  if (normalized_gene_filter$keep_low_or_high == "low") {
    keep_samples <- gene_expr < thr
  } else {
    keep_samples <- gene_expr > thr
  }

  sample_ids <- colnames(rnafilt_counts)[keep_samples]
  sample_ids <- trimws(as.character(sample_ids))
  unique(sample_ids[!is.na(sample_ids) & nzchar(sample_ids)])
}


matching_group_sample_ids <- function(clinic_annot,
                                      rnafilt_counts,
                                      clinic_filters = NULL,
                                      group_gene_filter = NULL,
                                      sample_id_col = "ID_Patient") {
  clinic_ids <- unique(extract_clinic_sample_ids(clinic_annot, sample_id_col = sample_id_col))

  if (!is.null(normalize_clinic_filters(clinic_filters))) {
    clinic_ids <- intersect(clinic_ids, matching_clinic_sample_ids(clinic_annot, clinic_filters, sample_id_col))
  }

  normalized_gene_filters <- normalize_group_gene_filters(group_gene_filter)
  if (!is.null(normalized_gene_filters)) {
    for (gene_filter in normalized_gene_filters) {
      clinic_ids <- intersect(clinic_ids, matching_gene_sample_ids(rnafilt_counts, gene_filter))
    }
  }

  clinic_ids
}


group_definition_path_suffix <- function(clinic_filters = NULL, group_gene_filter = NULL) {
  suffix_parts <- character(0)

  clinic_suffix <- clinic_filters_path_suffix(clinic_filters)
  if (!identical(clinic_suffix, "nofilt")) {
    suffix_parts <- c(suffix_parts, clinic_suffix)
  }

  gene_suffix <- group_gene_filters_path_suffix(group_gene_filter)
  if (!is.null(gene_suffix)) {
    suffix_parts <- c(suffix_parts, gene_suffix)
  }

  if (length(suffix_parts) == 0) {
    return("nofilt")
  }

  compact_path_component(paste(suffix_parts, collapse = "__"), max_chars = 120, prefix_chars = 88)
}


comparison_filters_path_suffix <- function(control_filters, test_filters, control_gene_filter = NULL, test_gene_filter = NULL) {
  compact_path_component(paste0(
    "control__", group_definition_path_suffix(control_filters, control_gene_filter),
    "__VS__",
    "test__", group_definition_path_suffix(test_filters, test_gene_filter)
  ), max_chars = 150, prefix_chars = 108)
}


filtering_on_clinic_and_genes <- function(clinic_annot = NULL, 
                                          rnafilt_counts = NULL, 
                                          output_dir = NULL,
                                          clinic_filters = NULL,
                                          filter_by_gene = NULL, 
                                          keep_low_or_high = NULL, 
                                          quantile_thr = NULL, 
                                          parsed_design = NULL,
                                          folder_name = "undefined_dir_name") 
{
  if (is.null(clinic_annot) & is.null(rnafilt_counts)){
    stop("No data to perform filtering on")
    return()
  }
  
  if (is.data.frame(rnafilt_counts)) {
    rnafilt_counts <- as.matrix(rnafilt_counts)
  }

  normalized_clinic_filters <- normalize_clinic_filters(clinic_filters)
  clinic_suffix <- clinic_filters_path_suffix(normalized_clinic_filters)
  parsed_design <- compact_path_component(parsed_design, max_chars = 170, prefix_chars = 120)

  # Filter using clinic annotations
  if (!is.null(normalized_clinic_filters)) {
    keep_samples <- rep(TRUE, nrow(clinic_annot))

    for (filter_name in names(normalized_clinic_filters)) {
      if (!filter_name %in% colnames(clinic_annot)) {
        stop(sprintf("Clinical column '%s' not found in clinic_annot.", filter_name))
      }

      keep_samples <- keep_samples & clinic_annot[, filter_name] %in% normalized_clinic_filters[[filter_name]]
    }

    clinic_annot <- clinic_annot[keep_samples, , drop = FALSE]
    filt_on_clin <- TRUE
  } else {
    filt_on_clin <- FALSE
  }

  # Filter using a gene expression level
  if (!is.null(filter_by_gene) &&
      filter_by_gene != "None" &&
      !is.null(quantile_thr) &&
      quantile_thr != 0 &&
      !is.null(keep_low_or_high) &&
      keep_low_or_high != "None") {

    if (!filter_by_gene %in% rownames(rnafilt_counts)) {
      stop(sprintf("Gene '%s' not found in rnafilt_counts.", filter_by_gene))
    }

    gene_expr <- as.numeric(rnafilt_counts[filter_by_gene, ])
    thr <- quantile(gene_expr, probs = quantile_thr, na.rm = TRUE)

    if (keep_low_or_high == "low") {
      keep_samples <- gene_expr < thr
    } else if (keep_low_or_high == "high") {
      keep_samples <- gene_expr > thr
    } else {
      stop("keep_low_or_high must be 'low' or 'high'")
    }

    cat("Samples kept after gene filter:", sum(keep_samples), "/", length(keep_samples), "\n")

    if (sum(keep_samples) == 0) {
      stop("Gene-based filtering removed all samples.")
    }

    rnafilt_counts <- rnafilt_counts[, keep_samples, drop = FALSE]
    filt_on_gene <- TRUE
  } else {
    filt_on_gene <- FALSE
  }

  # Elaborate the output path 
  gene_suffix <- NULL

  if (filt_on_gene) {
    gene_suffix <- compact_path_component(paste(
      sanitize_filter_token(filter_by_gene),
      sanitize_filter_token(keep_low_or_high),
      sanitize_filter_token(quantile_thr),
      sep = "-"
    ), max_chars = 80, prefix_chars = 60)
  }

  suffix_parts <- clinic_suffix
  if (!is.null(gene_suffix)) {
    suffix_parts <- compact_path_component(paste(c(suffix_parts, gene_suffix), collapse = "__"), max_chars = 120, prefix_chars = 88)
  }

  if (filt_on_gene && filt_on_clin) {
    output_path <- paste0(output_dir, "/", folder_name, "/", suffix_parts, "/", parsed_design)
  } else if (filt_on_clin) {
    output_path <- paste0(output_dir, "/", folder_name, "/", suffix_parts, "/", parsed_design)
  } else if (filt_on_gene) {
    output_path <- paste0(output_dir, "/", folder_name, "/", suffix_parts, "/", parsed_design)
  } else {
    output_path <- paste0(output_dir, "/", folder_name, "/", suffix_parts, "/", parsed_design)
  }
  
  return(list(
    clinic_annot = clinic_annot, 
    rnafilt_counts = rnafilt_counts, 
    output_path = output_path,
    clinic_filter_suffix = clinic_suffix
  ))
}
