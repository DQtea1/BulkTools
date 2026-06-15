library(jsonlite)
library(data.table)
library(dplyr)
### Contient ###
# normVST_bulk()
# prepare_projection()
# flatten_list_cols()
# symbol_to_ensg()


normVST_bulk <- function(rnaseq) {
    library(DESeq2)
    rnaseq <- as.matrix(rnaseq)

    # All-zero genes give a NaN geometric mean and add no information; drop them.
    rnaseq <- rnaseq[rowSums(rnaseq) > 0, , drop = FALSE]

    # An all-zero sample produces an NA size factor and cannot be transformed;
    # surface it clearly instead of failing deep inside DESeq2's sizeFactors<-.
    zero_samples <- colSums(rnaseq) == 0
    if (any(zero_samples)) {
        stop(sprintf(
            "normVST_bulk: %d sample(s) have zero counts across all genes and cannot be VST-normalized: %s",
            sum(zero_samples), paste(colnames(rnaseq)[zero_samples], collapse = ", ")
        ))
    }

    dds <- DESeqDataSetFromMatrix(countData = rnaseq,
                            colData = data.frame(cond=rnorm(ncol(rnaseq))>0),
                            design = ~ cond)

    # Default "ratio" size factors need at least one gene that is non-zero in
    # every sample. On sparse / unfiltered / union-merged matrices that may not
    # hold ("every gene contains at least one zero"), so fall back to the
    # zero-tolerant "poscounts" estimator. Dense matrices keep the default
    # behaviour (identical results); only the previously-failing case changes.
    dds <- tryCatch(
        estimateSizeFactors(dds),
        error = function(e) estimateSizeFactors(dds, type = "poscounts")
    )

    dds <- vst(dds, blind = TRUE)
    norm_matrix <- assay(dds)
    return(norm_matrix)
}


prepare_projection <- function(ref_bulk, sample_to_proj_filepath){
  library(tximport)
  samples <- data.frame(ID_Patient = "Studied_sample",
                        filepath = sample_to_proj_filepath,
                        stringsAsFactors = FALSE)
                        
  bulk_file <- samples$filepath
  names(bulk_file) <- samples$ID_Patient
  exemple_bulk = read.delim(bulk_file[1], sep = "\t")

  processed_sample = tximport(bulk_file, 
                          type = "salmon", 
                          txIn = TRUE, 
                          countsFromAbundance = c("lengthScaledTPM"), 
                          tx2gene = exemple_bulk[,1:2], 
                          geneIdCol = "Gene_name", 
                          txIdCol = "Name", 
                          abundanceCol = "TPM", 
                          countsCol = "NumReads", 
                          lengthCol = "Length", 
                          ignoreTxVersion = FALSE)

  studied_sample_counts <- processed_sample$counts
  ID_ref_samples <- colnames(ref_bulk)

  common_genes <- intersect(rownames(ref_bulk), rownames(studied_sample_counts))
  if (length(common_genes) == 0) {
    stop("Aucun gène en commun, vérifier les ID (ENSG, HUGO, .....)")
  }

  ref_sub <- ref_bulk[common_genes, , drop = FALSE]
  proj_vec <- studied_sample_counts[common_genes, 1, drop = TRUE]

  merged_bulks_counts <- cbind(ref_sub, Studied_sample = proj_vec)

  merged_bulks_vst <- normVST_bulk(round(merged_bulks_counts))

  return(list(merged_bulks_vst = merged_bulks_vst, 
              ID_ref_samples = ID_ref_samples, 
              sample_name = names(bulk_file)))
}


flatten_list_cols <- function(res, sep = ";", na_value = NA_character_) {
  DT <- as.data.table(res)

  # find list-type columns
  list_cols <- names(DT)[vapply(DT, is.list, logical(1))]

  # convert each list cell to a single string
  for (nm in list_cols) {
    DT[, (nm) := vapply(DT[[nm]], function(x) {
      if (is.null(x) || length(x) == 0) return(na_value)
      if (!is.atomic(x)) x <- unlist(x, recursive = TRUE, use.names = FALSE)
      paste(x, collapse = sep)
    }, character(1))]
  }

  DT
}




symbol_to_ensg <- function(mat, keep_unmapped = FALSE, collapse_duplicates = TRUE,
                          agg_fun = sum) {
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  stopifnot(!is.null(rownames(mat)))

  sym <- rownames(mat)

  # mapping SYMBOL -> ENSEMBL (ENSG...)
  ensg <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = sym,
    column = "ENSEMBL",
    keytype = "SYMBOL",
    multiVals = "first"   # si 1 symbole -> plusieurs ENSG, on prend le 1er
  )

  # Optionnel: essayer aussi les alias si beaucoup de NA
  # (ex: anciens symboles)
  if (anyNA(ensg)) {
    idx <- which(is.na(ensg))
    ensg2 <- AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = sym[idx],
      column = "ENSEMBL",
      keytype = "ALIAS",
      multiVals = "first"
    )
    ensg[idx] <- ensg2
  }

  if (!keep_unmapped) {
    keep <- !is.na(ensg)
    mat <- mat[keep, , drop = FALSE]
    ensg <- ensg[keep]
  } else {
    ensg[is.na(ensg)] <- sym[is.na(ensg)]  # on garde le symbole si pas mappé
  }

  rownames(mat) <- unname(ensg)

  # gérer les duplicats (plusieurs symboles -> même ENSG)
  if (collapse_duplicates) {
    # rowsum marche bien pour des counts; agg_fun = sum par défaut
    if (!identical(agg_fun, sum)) {
      # fallback général si tu veux mean/median etc.
      df <- data.frame(ENSG = rownames(mat), mat, check.names = FALSE)
      mat2 <- aggregate(. ~ ENSG, data = df, FUN = agg_fun)
      rn <- mat2$ENSG
      mat2$ENSG <- NULL
      mat2 <- as.matrix(mat2)
      rownames(mat2) <- rn
      return(mat2)
    } else {
      return(rowsum(mat, group = rownames(mat), reorder = FALSE))
    }
  }

  mat
}
