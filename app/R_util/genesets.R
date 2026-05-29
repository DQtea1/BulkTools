library(jsonlite)

## Load gene sets for ssGSEA, GSEA and/or GESECA ##
app_dir = normalizePath(getwd())
pathways_dir = file.path(app_dir, "REF_DATA", "genesets_and_pathways")

to_list = function(df) split(df$gene_symbol, df$gs_name)

normalize_gene_set_collection <- function(x, fallback_prefix = "geneset") {
  flatten_gene_sets <- function(obj, prefix = fallback_prefix) {
    if (!is.list(obj)) {
      genes <- trimws(as.character(obj))
      genes <- unique(genes[!is.na(genes) & nzchar(genes)])

      out <- list(genes)
      names(out) <- prefix
      return(out)
    }

    obj_names <- names(obj)
    if (is.null(obj_names)) {
      obj_names <- rep("", length(obj))
    }

    out <- list()

    for (idx in seq_along(obj)) {
      element_name <- trimws(as.character(obj_names[[idx]]))
      if (!nzchar(element_name)) {
        element_name <- paste0(prefix, "_", idx)
      }

      out <- c(out, flatten_gene_sets(obj[[idx]], prefix = element_name))
    }

    out
  }

  normalized <- flatten_gene_sets(x, prefix = fallback_prefix)
  keep <- vapply(normalized, length, integer(1)) > 0
  normalized <- normalized[keep]

  if (length(normalized) == 0) {
    return(list())
  }

  names(normalized) <- make.unique(names(normalized))
  normalized
}


all_pathways_raw = readRDS(file.path(pathways_dir, "human_pathways.rds"))
all_pathways = normalize_gene_set_collection(all_pathways_raw, fallback_prefix = "all_pathways")
c2_raw = fromJSON(file.path(pathways_dir, "c2_all_genesets.json"))

boyault_data = c2_raw[grepl("BOYAULT", names(c2_raw))]
boyault_sets = normalize_gene_set_collection(
  lapply(boyault_data, function(entry) entry[["geneSymbols"]]),
  fallback_prefix = "boyault"
)

c2_gene_sets = normalize_gene_set_collection(
  lapply(c2_raw, function(entry) entry[["geneSymbols"]]),
  fallback_prefix = "c2"
)


REACTOME_pathways = c2_gene_sets[grepl("REACTOME", names(c2_gene_sets), fixed = TRUE)]
KEGG_pathways = c2_gene_sets[grepl("KEGG", names(c2_gene_sets), fixed = TRUE)]
BIOCARTA_pathways = c2_gene_sets[grepl("BIOCARTA", names(c2_gene_sets), fixed = TRUE)]

GOBP_pathways = normalize_gene_set_collection(
  all_pathways_raw$c5[grepl("GOBP", names(all_pathways_raw$c5), fixed = TRUE)],
  fallback_prefix = "gobp"
)


transcript_factor_target = normalize_gene_set_collection(
  all_pathways_raw$c3[grepl("TARGET_GENES", names(all_pathways_raw$c3), fixed = TRUE)],
  fallback_prefix = "tft"
)

stress_msig_names = c(
  # proximal adrenergic / second messenger
  "GOBP_ADRENERGIC_RECEPTOR_SIGNALING_PATHWAY", "BIOCARTA_CREB_PATHWAY", "GOBP_RESPONSE_TO_MONOAMINE",
  "REACTOME_NOREPINEPHRINE_NEUROTRANSMITTER_RELEASE_CYCLE",
  # inflammatory / CTRA-adjacent
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  # stress-linked tumor outputs
  "HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_ANGIOGENESIS", "PID_LYMPH_ANGIOGENESIS_PATHWAY",
  # neural remodeling
  "GOBP_NEUROTROPHIN_TRK_RECEPTOR_SIGNALING_PATHWAY", "GOBP_AXON_GUIDANCE",
  # immune polarization / dysfunction
  # these DN sets are intentionally chosen because enrichment means "more M2" or "more exhausted"
  "COATES_MACROPHAGE_M1_VS_M2_DN", "GSE5099_CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_DN", "GSE9650_NAIVE_VS_EXHAUSTED_CD8_TCELL_DN", "GSE24026_PD1_LIGATION_VS_CTRL_IN_ACT_TCELL_LINE_UP"
)

paths_hallmark = intersect(stress_msig_names, names(all_pathways_raw$h))
paths_c2 = intersect(stress_msig_names, names(c2_gene_sets))
paths_c5 = intersect(stress_msig_names, names(all_pathways_raw$c5))

QoL_pathways = normalize_gene_set_collection(
  c(c2_gene_sets[paths_c2], all_pathways_raw$c5[paths_c5], all_pathways_raw$h[paths_hallmark]),
  fallback_prefix = "qol"
)
