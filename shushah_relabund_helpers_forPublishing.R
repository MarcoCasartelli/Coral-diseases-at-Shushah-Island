#### Shushah disease microbiome - shared relative-abundance helpers ####

# Sourced by all three project scripts so their phyloseq construction and
# relative-abundance maths can never drift apart 

# Steps
#   * phylogenetic tree attached to the phyloseq object BEFORE rarefaction (midpoint-rooted)
#   * clean filter: Kingdom == "Bacteria"; drop Family "Mitochondria" / Order "Chloroplast"
#   * drop global singletons (taxa_sums <= 1).
#   * rarefy to 42800, seed 69 (all 65 samples retained - min depth is 42850).
#   * per-sample relative abundance = taxon count / that sample's total count
#     (after agglomeration to the reporting rank).
#   * Family  -> plain tax_glom(taxrank = "Family"); every ASV unclassified at
#     family is summed into ONE "Unclassified" family.
#   * Genus   -> ASVs first tagged with a disambiguated "Taxon" label, then
#     agglomerated on that label
#   * group mean / min / max are taken over EVERY sample in the group (genus x
#     health-status), including samples where the taxon was not detected (0).

SHUSHAH_RAREFY_DEPTH <- 42800L
SHUSHAH_RAREFY_SEED  <- 69L
SHUSHAH_TAX_RANKS    <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

# --- phyloseq construction -------------------------------------------------

# Build a phyloseq object from the four standard project files
shushah_build_phyloseq <- function(counts_file, taxonomy_file, metadata_file,
                                   tree_file, metadata_prep = NULL) {
  taxonomy <- as.matrix(read.csv(taxonomy_file, header = TRUE, row.names = 1, sep = "\t"))
  counts   <- as.matrix(read.csv(counts_file,   header = TRUE, row.names = 1, sep = "\t"))
  metadata <- read.csv(metadata_file, header = TRUE, row.names = 1)
  metadata <- data.frame(metadata, check.names = FALSE)
  if (!"SampleID" %in% colnames(metadata)) metadata$SampleID <- rownames(metadata)
  if (!is.null(metadata_prep)) {
    metadata <- metadata_prep(metadata)
    stopifnot(is.data.frame(metadata))
  }

  tree <- phyloseq::read_tree(tree_file)
  tree <- phytools::midpoint.root(tree)

  phyloseq::phyloseq(
    phyloseq::otu_table(counts, taxa_are_rows = TRUE),
    phyloseq::tax_table(taxonomy),
    phyloseq::sample_data(metadata),
    phyloseq::phy_tree(tree)
  )
}

# Clean filter, identical across all three scripts.
shushah_clean_phyloseq <- function(ps) {
  ps <- phyloseq::subset_taxa(
    ps, Kingdom == "Bacteria" & Family != "Mitochondria" & Order != "Chloroplast"
  )
  phyloseq::prune_taxa(phyloseq::taxa_sums(ps) > 1, ps)
}

# Rarefy to the fixed project depth / seed. set.seed() as well as rngseed=
shushah_rarefy <- function(ps, depth = SHUSHAH_RAREFY_DEPTH, seed = SHUSHAH_RAREFY_SEED,
                           verbose = FALSE) {
  set.seed(seed)
  phyloseq::rarefy_even_depth(
    ps, sample.size = depth, replace = FALSE, rngseed = seed, verbose = verbose
  )
}

# --- disambiguated genus tagging ----------------------------------------------

# Add a "Taxon" column to tax_table(ps): the ASV's own value at `level` if
# classified, otherwise its full upstream lineage pasted with "_Unclassified".
shushah_add_taxon_column <- function(ps, level = "Genus") {
  tt <- as.data.frame(as(phyloseq::tax_table(ps), "matrix"), stringsAsFactors = FALSE)
  rank_idx <- which(SHUSHAH_TAX_RANKS == level)
  if (length(rank_idx) == 0) {
    stop("shushah_add_taxon_column(): level must be one of ",
         paste(SHUSHAH_TAX_RANKS, collapse = ", "), call. = FALSE)
  }
  upstream <- SHUSHAH_TAX_RANKS[seq_len(rank_idx - 1)]
  lineage  <- if (length(upstream) > 0) {
    do.call(paste, c(tt[upstream], sep = "_"))
  } else {
    rep("Root", nrow(tt))
  }
  unclassified <- is.na(tt[[level]]) | tt[[level]] %in% c("Unclassified", "")
  tt$Taxon <- ifelse(unclassified, paste(lineage, "Unclassified", sep = "_"), tt[[level]])

  new_tax <- as.matrix(tt)
  rownames(new_tax) <- phyloseq::taxa_names(ps)
  phyloseq::tax_table(ps) <- phyloseq::tax_table(new_tax)
  ps
}

# --- per-sample relative abundance (long) -----------------------------------

.shushah_consensus_ranks <- function(tt, key, upstream_cols, taxon_col) {
  parts <- split(tt[, upstream_cols, drop = FALSE], key)
  rows  <- lapply(parts, function(g) {
    vapply(g, function(col) {
      u <- unique(as.character(col))
      if (length(u) == 1L) u else NA_character_
    }, character(1))
  })
  out <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  out[[taxon_col]] <- names(rows)
  tibble::as_tibble(out[, c(taxon_col, upstream_cols), drop = FALSE])
}

# Long per-sample relative-abundance table for one reporting rank.
#   level = "Family" -> key = the Family string; every ASV unclassified at
#                       family is summed into ONE "Unclassified" family.
#   level = "Genus"  -> key = disambiguated "Taxon" label (own genus name if
#                       classified, else full upstream lineage + "_Unclassified").
# Returns a tibble: <key col> + upstream rank columns + SampleID +
# RelativeAbundance (proportion 0-1; denominator = that sample's total count).
shushah_relabund_long <- function(ps, level = c("Family", "Genus")) {
  level <- match.arg(level)

  counts <- as(phyloseq::otu_table(ps), "matrix")
  if (!phyloseq::taxa_are_rows(ps)) counts <- t(counts)          # taxa x samples

  if (level == "Family") {
    ps_tag    <- ps
    taxon_col <- "Family"
    upstream  <- c("Kingdom", "Phylum", "Class", "Order")
  } else {
    ps_tag    <- shushah_add_taxon_column(ps, level = "Genus")
    taxon_col <- "Taxon"
    upstream  <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
  }
  tt <- as.data.frame(as(phyloseq::tax_table(ps_tag), "matrix"), stringsAsFactors = FALSE)
  key <- tt[[taxon_col]]
  key[is.na(key) | key == ""] <- "Unclassified"

  agg  <- rowsum(counts, group = key)                            # unique keys x samples
  prop <- sweep(agg, 2, colSums(agg), "/")

  rank_map <- .shushah_consensus_ranks(tt, key, upstream, taxon_col)

  long <- tibble::as_tibble(as.data.frame(prop), rownames = taxon_col)
  long <- tidyr::pivot_longer(long, -dplyr::all_of(taxon_col),
                              names_to = "SampleID", values_to = "RelativeAbundance")
  dplyr::left_join(rank_map, long, by = taxon_col)
}

# --- group mean / min / max ------------------------------------------------

# Collapse a shushah_relabund_long() table to one row per taxon per group.
#   long_df      : output of shushah_relabund_long()
#   sample_meta  : data.frame / tibble with SampleID + the grouping columns
#   group_cols   : character vector, e.g. c("GENUS", "STATUS_ABBR")
# mean/min/max are over EVERY sample in the group (zeros for non-detections included)
shushah_group_summary <- function(long_df, sample_meta, group_cols) {
  taxon_col <- if ("Taxon" %in% names(long_df)) "Taxon" else "Family"
  rank_cols <- intersect(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Taxon"),
                         names(long_df))
  rank_map  <- dplyr::distinct(long_df[, unique(c(taxon_col, rank_cols)), drop = FALSE])

  sample_meta <- tibble::as_tibble(sample_meta)[, unique(c("SampleID", group_cols)), drop = FALSE]

  df <- dplyr::left_join(long_df, sample_meta, by = "SampleID")
  df <- df[stats::complete.cases(df[, group_cols, drop = FALSE]), , drop = FALSE]

  summ <- dplyr::summarise(
    dplyr::group_by(df, dplyr::across(dplyr::all_of(c(group_cols, taxon_col)))),
    n_samples              = dplyr::n(),
    n_present              = sum(RelativeAbundance > 0),
    mean_RelativeAbundance = mean(RelativeAbundance),
    min_RelativeAbundance  = min(RelativeAbundance),
    max_RelativeAbundance  = max(RelativeAbundance),
    .groups = "drop"
  )
  summ <- dplyr::left_join(summ, rank_map, by = taxon_col)
  dplyr::arrange(summ, dplyr::across(dplyr::all_of(group_cols)),
                 dplyr::desc(mean_RelativeAbundance))
}

# --- top-N (per subset) UNION + "Others" ----------------------------------

# Collapse a shushah_relabund_long() table to the "top-N + Others" view used in
# the manuscript composition table.
shushah_topN_union_summary <- function(long_df, sample_meta, group_cols,
                                       subset_col, top_n = 19,
                                       others_label = "Others") {
  taxon_col <- if ("Taxon" %in% names(long_df)) "Taxon" else "Family"
  rank_cols <- intersect(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Taxon"),
                         names(long_df))
  rank_map  <- dplyr::distinct(long_df[, unique(c(taxon_col, rank_cols)), drop = FALSE])

  sample_meta <- tibble::as_tibble(sample_meta)[, unique(c("SampleID", group_cols, subset_col)), drop = FALSE]

  df <- dplyr::left_join(long_df, sample_meta, by = "SampleID")
  df <- df[stats::complete.cases(df[, unique(c(group_cols, subset_col)), drop = FALSE]), , drop = FALSE]

  # 1-2. top-N per subset, by mean per-sample rel. abundance; union of names
  top_by_subset <- df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(subset_col, taxon_col)))) %>%
    dplyr::summarise(.mean = mean(RelativeAbundance), .groups = "drop_last") %>%
    dplyr::slice_max(.mean, n = top_n, with_ties = FALSE) %>%
    dplyr::ungroup()
  named <- sort(unique(top_by_subset[[taxon_col]]))

  # 3. relabel non-union taxa as "Others", re-aggregate per (sample x label)
  df$.label <- ifelse(df[[taxon_col]] %in% named, df[[taxon_col]], others_label)
  per_sample <- df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c("SampleID", group_cols, ".label")))) %>%
    dplyr::summarise(RelativeAbundance = sum(RelativeAbundance), .groups = "drop")
  # make every (sample x label) cell explicit so absent named taxa contribute 0
  per_sample <- per_sample %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c("SampleID", group_cols)))) %>%
    tidyr::complete(.label = c(named, others_label),
                    fill = list(RelativeAbundance = 0)) %>%
    dplyr::ungroup()

  # 4. per-group mean / min / max
  summ <- per_sample %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(group_cols, ".label")))) %>%
    dplyr::summarise(
      n_samples              = dplyr::n(),
      mean_RelativeAbundance = mean(RelativeAbundance),
      min_RelativeAbundance  = min(RelativeAbundance),
      max_RelativeAbundance  = max(RelativeAbundance),
      .groups = "drop"
    )
  names(summ)[names(summ) == ".label"] <- taxon_col

  # order rows: named taxa by overall mean (desc), "Others" always last
  taxon_order <- summ %>%
    dplyr::filter(.data[[taxon_col]] != others_label) %>%
    dplyr::group_by(.data[[taxon_col]]) %>%
    dplyr::summarise(.m = mean(mean_RelativeAbundance), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(.m)) %>% dplyr::pull(1)
  summ[[taxon_col]] <- factor(summ[[taxon_col]], levels = c(taxon_order, others_label))
  summ <- dplyr::arrange(summ, dplyr::across(dplyr::all_of(group_cols)), summ[[taxon_col]])
  summ[[taxon_col]] <- as.character(summ[[taxon_col]])

  summ <- dplyr::left_join(summ, rank_map, by = taxon_col)   # ranks NA for "Others"
  attr(summ, "named_taxa") <- named
  summ
}

# --- methods note --------------------------------------------------------------

# Write the RELAB_METHODS.txt companion file next to an exported table set.
shushah_write_relab_methods <- function(path, use_rarefied, group_desc,
                                        levels_desc = "Family and Genus",
                                        extra_lines = NULL) {
  lines <- c(
    "Relative abundance - computation method (Shushah disease microbiome)",
    "===================================================================",
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("Helper: ShushahDisease_RelAb_Helpers/shushah_relabund_helpers.R"),
    "",
    "Phylogenetic tree attached to the phyloseq object BEFORE rarefaction",
    "(midpoint-rooted). Cleaning filter, applied in this order:",
    "  1. keep Kingdom == 'Bacteria'",
    "  2. drop Family == 'Mitochondria'; drop Order == 'Chloroplast'",
    "  3. drop global singletons (taxa with total count <= 1 across all samples)",
    "",
    if (isTRUE(use_rarefied)) {
      sprintf(paste0("Rarefaction: rarefy_even_depth(sample.size = %d, seed = %d, ",
                     "replace = FALSE). All 65 samples retained (minimum depth 42850)."),
              SHUSHAH_RAREFY_DEPTH, SHUSHAH_RAREFY_SEED)
    } else {
      "Rarefaction: NOT applied - relative abundances computed on the cleaned, non-rarefied counts."
    },
    "",
    "Per-sample relative abundance = taxon count / that sample's total count,",
    "computed after agglomeration to the reporting rank.",
    "",
    "Agglomeration (ASV counts summed onto the reporting key, then divided by",
    "the per-sample total):",
    "  Family - key = the Family string; every ASV unclassified at family is",
    "           summed into a single 'Unclassified' family.",
    "  Genus  - key = a disambiguated 'Taxon' label (own genus name if",
    "           classified, else full upstream lineage + '_Unclassified').",
    "",
    paste0("Group mean / min / max grouping: ", group_desc),
    "mean / min / max are taken over EVERY sample in the group, including samples",
    "where the taxon was not detected (counted as 0).",
    "",
    "Files per level:",
    "  <level>_relabundance_persample.csv       - one row per taxon x sample.",
    "  <level>_relabundance_group_summary.csv   - one row per taxon per group",
    "                                            (all taxa, no collapsing).",
    "  <level>_relabundance_top19_group_summary.csv - collapsed manuscript view:",
    "     top 19 taxa are picked WITHIN each coral genus (Acropora, Pocillopora)",
    "     by mean per-sample relative abundance across all that genus's samples;",
    "     the named set is the UNION of the two genera's top-19 lists; every",
    "     other taxon is folded into one 'Others' row (so 'Others' excludes the",
    "     top-19 of BOTH genera). Per group, the mean_RelativeAbundance of all",
    "     rows (named taxa + Others) sums to exactly 1 (100%); min / max are",
    "     per-row and do NOT sum to 1.",
    "",
    paste0("Levels exported: ", levels_desc),
    "Units: RelativeAbundance is a proportion in [0, 1] (multiply by 100 for %)."
  )
  if (!is.null(extra_lines)) lines <- c(lines, "", extra_lines)
  writeLines(lines, path)
  invisible(path)
}
