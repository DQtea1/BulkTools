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


clinic_filters_path_suffix <- function(clinic_filters) {
  normalized_filters <- normalize_clinic_filters(clinic_filters)

  if (is.null(normalized_filters)) {
    return("nofilt")
  }

  suffix_parts <- vapply(names(normalized_filters), function(filter_name) {
    values <- sort(unique(normalized_filters[[filter_name]]))
    paste0(
      sanitize_filter_token(filter_name),
      "-",
      paste(sanitize_filter_token(values), collapse = "+")
    )
  }, character(1))

  paste(suffix_parts, collapse = "__")
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

  message("normalized_clinic_filters: ", normalized_clinic_filters)
  message("clinic_suffix: ", clinic_suffix)

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

  message("clinic_annot : ", clinic_annot )
  message("keep_samples : ", keep_samples)

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
    gene_suffix <- paste(
      sanitize_filter_token(filter_by_gene),
      sanitize_filter_token(keep_low_or_high),
      sanitize_filter_token(quantile_thr),
      sep = "-"
    )
  }

  suffix_parts <- clinic_suffix
  if (!is.null(gene_suffix)) {
    suffix_parts <- paste(c(suffix_parts, gene_suffix), collapse = "__")
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

  message("output_path: ", output_path)
  
  return(list(
    clinic_annot = clinic_annot, 
    rnafilt_counts = rnafilt_counts, 
    output_path = output_path,
    clinic_filter_suffix = clinic_suffix
  ))
}
