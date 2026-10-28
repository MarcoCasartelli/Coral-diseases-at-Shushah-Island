#### ANCOM-BC2 Differential Abundance Pipeline ####
#### Libraries ####

library(tidyverse)   # data wrangling / ggplot2
library(phyloseq)    # microbiome data object + manipulation
library(ANCOMBC)     # ancombc2()
library(patchwork)   # combining the final figure
library(openxlsx)    # approved family-color palette

# Shared relative-abundance helpers to compute relative abundances consistently across codes
source("./shushah_relabund_helpers.R") #Edit with correct path

#### Configuration ####

setwd(
  "./ShushahDisease_Microbiome_PresenceAbsence" #Edit with correct path
)
out_dir <- "ANCOMBC2_Results"
dir.create(out_dir, showWarnings = FALSE)

tax_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

## Fixed bubble-size from the presence-absence script, to keep sizes comparable
bubble_size_max <- 0.211

RELAB_USE_RAREFY <- TRUE
## TRUE  → Mean/Min/Max relative abundance values and plots are generated on the rarefied datasets
## FALSE → Uses cleaned, non rarefied dataset for the values above
## Note: ANCOM-BC2 is always run on the non-rarefied object, regardless of this switch

#### Import + build ASV-level phyloseq (no rarefaction - ANCOM-BC2 does its own bias correction) ####

taxfile <- as.matrix(read.csv("ShushahMicrobiome_Total_ASVs_taxonomy-no-contam.tsv", header = TRUE, row.names = 1, sep = "\t"))
mapfile <- read.csv("ShushahMicrobiome_Total_Metadata.csv", header = TRUE, row.names = 1)
ASVs <- as.matrix(read.csv("ShushahMicrobiome_Total_ASVs_counts-no-contam.tsv", header = TRUE, row.names = 1, sep = "\t"))

mapfile <- data.frame(mapfile)
mapfile$SampleID <- rownames(mapfile)
mapfile <- mapfile %>% unite(col = "STATUS", COLONY_STATUS, FRAGMENT_STATUS, sep = "_", remove = FALSE)
status_levels <- c("Healthy_NoLesion", "Diseased_NoLesion", "Diseased_Lesion")
status_abbr   <- c("AH", "DU", "DL")
mapfile$STATUS <- factor(mapfile$STATUS, levels = status_levels)
mapfile$STATUS_ABBR <- factor(
  dplyr::recode(mapfile$STATUS, "Healthy_NoLesion" = "AH", "Diseased_NoLesion" = "DU", "Diseased_Lesion" = "DL"),
  levels = status_abbr
)
mapfile$MONTH <- factor(mapfile$MONTH, levels = c("March", "July"))

## The tree is attached to maintain relative abundances comparable across scripts
tree_ANC <- phytools::midpoint.root(
  phyloseq::read_tree("ShushahMicrobiome_Total_ASVs-no-contam_tree.nwk")
)
ps <- phyloseq(otu_table(ASVs, taxa_are_rows = TRUE), tax_table(taxfile),
               sample_data(mapfile), phy_tree(tree_ANC))

## Same organelle/singleton filter as the other scripts' pipelines
## (shared helper - identical across all three scripts).
ps_clean <- shushah_clean_phyloseq(ps)
ps_clean

#### Disambiguated genus-level tagging ####
tag_taxon_column <- function(ps) {
  tt_df <- as.data.frame(as(tax_table(ps), "matrix"), stringsAsFactors = FALSE)
  upstream_ranks <- tax_ranks[seq_len(which(tax_ranks == "Genus") - 1)]
  lineage <- do.call(paste, c(tt_df[upstream_ranks], sep = "_"))
  tt_df$Taxon <- ifelse(!is.na(tt_df$Genus) & tt_df$Genus != "Unclassified", tt_df$Genus, paste(lineage, "Unclassified", sep = "_"))
  new_tax <- as.matrix(tt_df)
  rownames(new_tax) <- taxa_names(ps)
  tax_table(ps) <- tax_table(new_tax)
  ps
}
ps_clean <- tag_taxon_column(ps_clean)

## Number of ASVs that originally belong to each disambiguated genus ("Taxon"),
## counted across the FULL cleaned (non-rarefied, all 65 samples) dataset
taxon_asv_n <- as.data.frame(
  table(Taxon = as(tax_table(ps_clean), "matrix")[, "Taxon"]),
  responseName = "n_ASVs", stringsAsFactors = FALSE
)
taxon_asv_n$Taxon  <- as.character(taxon_asv_n$Taxon)
taxon_asv_n$n_ASVs <- as.integer(taxon_asv_n$n_ASVs)

## Readable display label for a Taxon value
make_display_label <- function(df) {
  ranks_deep_to_shallow <- rev(tax_ranks)
  label <- rep(NA_character_, nrow(df))
  matched_rank <- rep(NA_character_, nrow(df))
  for (r in ranks_deep_to_shallow) {
    todo <- is.na(label)
    hit <- todo & !is.na(df[[r]]) & df[[r]] != "Unclassified"
    label[hit] <- df[[r]][hit]
    matched_rank[hit] <- r
  }
  ifelse(matched_rank == "Genus", label, paste0(label, " (", tolower(matched_rank), ")"))
}

## Family-color palette
family_palette_file <- "./Color codes families SVG.xlsx" #Edit with correct path
approved_colors <- read.xlsx(family_palette_file, sheet = 1) %>%
  mutate(HexColor = paste0("#", substr(Color, 1, 6))) %>%
  { setNames(.$HexColor, .$Family) }

cat("Setup complete:", ntaxa(ps_clean), "taxa,", nsamples(ps_clean), "samples,",
    length(unique(as(tax_table(ps_clean), "matrix")[, "Taxon"])), "disambiguated genus-level groups.\n")

#### Rarefaction (for relative abundance and plot only) ####
## Same depth/seed as the rest of the project (shared helper)
rarefy_depth <- 42800
rarefy_seed <- 69
ps_rarefied <- shushah_rarefy(ps_clean, depth = rarefy_depth, seed = rarefy_seed)

## Relative-abundance source object (see RELAB_USE_RAREFY in Configuration).
ps_relab_src <- if (isTRUE(RELAB_USE_RAREFY)) ps_rarefied else ps_clean

## Genus-level (disambiguated Taxon) relative abundance
taxon_relabund <- function(ps) {
  taxon <- as(tax_table(ps), "matrix")[, "Taxon"]
  counts <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) counts <- t(counts)
  agg <- rowsum(counts, group = taxon)
  sweep(agg, 2, colSums(agg), "/")
}

#### Relative abundance exports ####
local({
  relab_dir <- file.path(out_dir, "relative_abundance")
  dir.create(relab_dir, showWarnings = FALSE)

  sm <- as(sample_data(ps_relab_src), "data.frame") %>%
    transmute(SampleID, GENUS = as.character(GENUS), STATUS_ABBR)

  for (lvl in c("Family", "Genus")) {
    long  <- shushah_relabund_long(ps_relab_src, level = lvl)
    summ  <- shushah_group_summary(long, sm, group_cols = c("GENUS", "STATUS_ABBR"))
    top19 <- shushah_topN_union_summary(long, sm, group_cols = c("GENUS", "STATUS_ABBR"),
                                        subset_col = "GENUS", top_n = 19)
    write.csv(long,  file.path(relab_dir, paste0(lvl, "_relabundance_persample.csv")),
              row.names = FALSE)
    write.csv(summ,  file.path(relab_dir, paste0(lvl, "_relabundance_group_summary.csv")),
              row.names = FALSE)
    write.csv(top19, file.path(relab_dir, paste0(lvl, "_relabundance_top19_group_summary.csv")),
              row.names = FALSE)
    cat("  ", lvl, ": ", nrow(long), " per-sample rows, ", nrow(summ),
        " group rows, ", length(attr(top19, "named_taxa")),
        " named + Others (top-19 union)\n", sep = "")
  }

  shushah_write_relab_methods(
    file.path(relab_dir, "RELAB_METHODS.txt"),
    use_rarefied = isTRUE(RELAB_USE_RAREFY),
    group_desc   = paste("coral genus x health-status (AH = Healthy/NoLesion,",
                         "DU = Diseased/NoLesion, DL = Diseased/Lesion); March + July pooled"),
    extra_lines  = c(
      "Source script: ShushahTotal_Microbiome_ANCOMBC2_forPublishing.R",
      "Relative-abundance summaries only. The ANCOM-BC2 differential-abundance",
      "test itself runs on the raw (non-rarefied) counts and is unaffected by",
      "RELAB_USE_RAREFY."
    )
  )
  cat("Relative abundance tables written to ", relab_dir, "/\n", sep = "")
})

#### Helper: run one within-genus pairwise ANCOM-BC2 comparison ####
## group_a/group_b: two of "Healthy_NoLesion"/"Diseased_NoLesion"/"Diseased_Lesion" 
## group_a is the reference level, so a positive lfc means enriched in group_b.
run_comparison <- function(ps_raw, ps_rar, coral_genus, group_a, group_b) {
  smp <- as(sample_data(ps_raw), "data.frame")
  keep <- rownames(smp)[smp$GENUS == coral_genus & smp$STATUS %in% c(group_a, group_b)]
  ps_sub <- prune_taxa(taxa_sums(prune_samples(keep, ps_raw)) > 0, prune_samples(keep, ps_raw))
  sample_data(ps_sub)$STATUS <- ordered(sample_data(ps_sub)$STATUS, levels = c(group_a, group_b))

  fit <- ancombc2(
    ps_sub, fix_formula = "STATUS", tax_level = "Taxon", p_adj_method = "BH",
    lib_cut = 1000, group = "STATUS", struc_zero = TRUE, neg_lb = FALSE, alpha = 0.05,
    em_control = list(tol = 1e-05, max_iter = 100)
  )
  res <- fit$res

  sig_robust <- res %>%
    filter(diff_robust_STATUS.L == TRUE) %>%
    transmute(
      Taxon = taxon, lfc = lfc_STATUS.L, se = se_STATUS.L,
      w_statistic = W_STATUS.L, q = q_STATUS.L,
      wg = ifelse(lfc_STATUS.L > 0, "pos", "neg")
    )

  ## Relative abundance (rarefied), restricted to the same sample subset
  keep_rar <- intersect(keep, sample_names(ps_rar))
  ps_sub_rar <- prune_taxa(taxa_sums(prune_samples(keep_rar, ps_rar)) > 0, prune_samples(keep_rar, ps_rar))
  relab <- taxon_relabund(ps_sub_rar)

  abbr_a <- status_abbr[match(group_a, status_levels)]
  abbr_b <- status_abbr[match(group_b, status_levels)]

  ## Results summary for each robust genus
  smp_rar  <- as(sample_data(ps_sub_rar), "data.frame")
  cols_by  <- list(); cols_by[[abbr_a]] <- rownames(smp_rar)[smp_rar$STATUS == group_a]
  cols_by[[abbr_b]] <- rownames(smp_rar)[smp_rar$STATUS == group_b]
  grp_stats <- function(tx) {
    out <- list()
    for (ab in c(abbr_a, abbr_b)) {
      cols   <- cols_by[[ab]]
      has    <- tx %in% rownames(relab)
      v      <- if (has) relab[tx, cols] else rep(NA_real_, length(cols))
      n_pres <- if (has) sum(v > 0) else NA_integer_
      out <- c(out, setNames(
        list(n_pres, length(cols),
             if (has) n_pres / length(cols) else NA_real_,
             mean(v), suppressWarnings(min(v)), suppressWarnings(max(v))),
        paste0(c("N_SamplesPresent_", "N_SamplesInGroup_", "Prevalence_",
                 "mean_RelAb_", "min_RelAb_", "max_RelAb_"), ab)))
    }
    as.data.frame(out, check.names = FALSE)
  }

  ## ASV counts per robust genus:
  ##   n_ASVs          - in the FULL cleaned dataset (all 65 samples, non-rarefied);
  ##   n_ASVs_contrast - ASVs with >=1 read ANYWHERE in this comparison's sample set (ps_sub)
  ##   n_ASVs_<group>  - of the genus's ASVs, how many are detected (>=1 read) in
  ##                     each of the two contrast groups' samples.
  otu_sub <- as(otu_table(ps_sub), "matrix"); if (!taxa_are_rows(ps_sub)) otu_sub <- t(otu_sub)
  tx_sub  <- as(tax_table(ps_sub), "matrix")[, "Taxon"]
  smp_sub <- as(sample_data(ps_sub), "data.frame")
  ca <- rownames(smp_sub)[smp_sub$STATUS == group_a]
  cb <- rownames(smp_sub)[smp_sub$STATUS == group_b]
  pres_a <- rowSums(otu_sub[, ca, drop = FALSE] > 0) > 0
  pres_b <- rowSums(otu_sub[, cb, drop = FALSE] > 0) > 0
  asv_contrast_df <- data.frame(
    Taxon           = names(table(tx_sub)),
    n_ASVs_contrast = as.integer(table(tx_sub)),
    stringsAsFactors = FALSE
  )
  asv_contrast_df[[paste0("n_ASVs_", abbr_a)]] <- as.integer(tapply(pres_a, tx_sub, sum)[asv_contrast_df$Taxon])
  asv_contrast_df[[paste0("n_ASVs_", abbr_b)]] <- as.integer(tapply(pres_b, tx_sub, sum)[asv_contrast_df$Taxon])

  if (nrow(sig_robust) > 0) {
    stats_df <- do.call(rbind, lapply(sig_robust$Taxon, grp_stats))
    per_group_block <- function(ab) {
      paste0(c("n_ASVs_", "N_SamplesPresent_", "N_SamplesInGroup_", "Prevalence_",
               "mean_RelAb_", "min_RelAb_", "max_RelAb_"), ab)
    }
    sig_robust <- dplyr::bind_cols(sig_robust, stats_df) %>%
      dplyr::left_join(taxon_asv_n,     by = "Taxon") %>%
      dplyr::left_join(asv_contrast_df, by = "Taxon") %>%
      dplyr::relocate(dplyr::all_of(c("n_ASVs", "n_ASVs_contrast",
                                      per_group_block(abbr_a),
                                      per_group_block(abbr_b))),
                      .after = wg)
  }

  taxonomy <- as.data.frame(as(tax_table(ps_sub), "matrix"), stringsAsFactors = FALSE) %>%
    distinct(Taxon, .keep_all = TRUE) %>%
    mutate(TaxonLabel = make_display_label(.))

  list(
    coral_genus = coral_genus, group_a = group_a, group_b = group_b,
    n_tested = nrow(res), n_sig = sum(res$diff_STATUS.L, na.rm = TRUE), n_robust = nrow(sig_robust),
    sig_robust = sig_robust, relab = relab, taxonomy = taxonomy,
    sample_ids = colnames(relab)
  )
}

#### Comparison's plotting data (bar + bubble) ####

## comparison_title: short text for axis/plot titles, e.g. "AH vs DL".
build_comparison_data <- function(cmp, comparison_title) {
  bar_df <- cmp$sig_robust %>%
    left_join(cmp$taxonomy, by = "Taxon") %>%
    mutate(TaxonLabel = make.unique(TaxonLabel, sep = " #"))

  ## Row order, TOP to BOTTOM: descending lfc 
  taxon_order <- bar_df %>% arrange(desc(lfc)) %>% pull(TaxonLabel) %>% rev()
  bar_df <- bar_df %>% mutate(TaxonLabel = factor(TaxonLabel, levels = taxon_order))

  ## Family palette
  bar_df <- bar_df %>% mutate(FamilyColorGroup = ifelse(Family %in% names(approved_colors), Family, "Others"))

  ## Bubble data: one row per (taxon, sample)
  relab <- cmp$relab[as.character(bar_df$Taxon), , drop = FALSE]
  bubble_df <- as_tibble(relab, rownames = "Taxon") %>%
    pivot_longer(-Taxon, names_to = "SampleID", values_to = "RelativeAbundance") %>%
    mutate(RelativeAbundance = ifelse(RelativeAbundance > 0, RelativeAbundance, NA)) %>%
    left_join(bar_df %>% select(Taxon, TaxonLabel, FamilyColorGroup), by = "Taxon") %>%
    left_join(as(sample_data(ps_rarefied), "data.frame") %>% select(SampleID, STATUS_ABBR, MONTH), by = "SampleID")

  ## Factor reorder
  sample_order <- as(sample_data(ps_rarefied), "data.frame") %>%
    filter(SampleID %in% cmp$sample_ids) %>%
    mutate(STATUS_ABBR = factor(STATUS_ABBR, levels = status_abbr), MONTH = factor(MONTH, levels = c("March", "July"))) %>%
    arrange(MONTH, STATUS_ABBR, SampleID) %>%
    mutate(SampleLabel = paste(SampleID, STATUS_ABBR, MONTH, sep = "\n"))
  sample_order$SampleLabel <- factor(sample_order$SampleLabel, levels = sample_order$SampleLabel)
  bubble_df <- bubble_df %>%
    left_join(sample_order %>% select(SampleID, SampleLabel), by = "SampleID")

  color_group_levels <- bar_df %>%
    mutate(ColorSortKey = ifelse(FamilyColorGroup == "Others", "zzz_Others", FamilyColorGroup)) %>%
    distinct(FamilyColorGroup, ColorSortKey) %>%
    arrange(ColorSortKey) %>%
    pull(FamilyColorGroup)

  list(
    title = comparison_title, bar_df = bar_df, bubble_df = bubble_df, sample_order = sample_order,
    taxon_order = taxon_order, n_taxa = length(taxon_order), n_samples = nrow(sample_order),
    color_group_levels = color_group_levels, family_palette = approved_colors[color_group_levels]
  )
}

#### Plot functions ####

## Log-fold-change bar plot 
## lfc_limits/bar_width are shared across every comparison in the figure
plot_lfc_bar <- function(cd, lfc_limits, bar_width) {
  ggplot(cd$bar_df, aes(x = TaxonLabel, y = lfc, fill = wg)) +
    geom_col(width = bar_width) +
    geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), width = bar_width * 0.4, linewidth = 0.3, color = "black") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    coord_flip(ylim = lfc_limits) +
    scale_fill_manual(values = c(pos = "grey30", neg = "grey75"), guide = "none") +
    scale_y_continuous(expand = expansion(mult = 0.02)) +
    labs(x = NULL, y = paste0("Log fold change (", cd$title, ")")) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.y = element_text(face = "italic", size = 8),
      axis.text.x = element_text(size = 8),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 11)
    )
}

## Bubble plot: 
## fill = Family
## size = relative abundance calibrated against the fixed shared reference bubble_size_max
plot_bubble <- function(cd) {
  cd$bubble_df <- cd$bubble_df %>% mutate(FamilyColorGroup = factor(FamilyColorGroup, levels = cd$color_group_levels))
  ggplot(cd$bubble_df, aes(x = SampleLabel, y = TaxonLabel, size = RelativeAbundance, fill = FamilyColorGroup)) +
    geom_point(shape = 21, color = "grey30", stroke = 0.15, alpha = 0.9) +
    scale_size_area(max_size = 8, limits = c(0, bubble_size_max), oob = function(x, range) x,
                     breaks = scales::extended_breaks(n = 5)(c(0, bubble_size_max)),
                     labels = scales::percent_format(accuracy = 0.1), name = "Relative\nabundance") +
    scale_fill_manual(values = cd$family_palette, name = "Family", drop = FALSE) +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    guides(fill = guide_legend(override.aes = list(size = 4), ncol = 1)) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major = element_line(color = "grey93"),
      panel.grid.minor = element_blank()
    )
}

## Thin month-shape strip: black triangle (March) / black circle (July),
plot_month_strip <- function(cd) {
  ggplot(cd$sample_order, aes(x = SampleLabel, y = 1, shape = MONTH)) +
    geom_point(size = 1.8, color = "black") +
    scale_shape_manual(values = c(March = 17, July = 16), guide = "none") +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(limits = c(0.9, 1.1)) +
    labs(x = NULL, y = NULL) +
    theme_void(base_size = 10) +
    theme(
      axis.text.x = element_text(size = 6, lineheight = 0.8, color = "black", margin = margin(t = 2)),
      plot.margin = margin(0, 5.5, 0, 5.5)
    )
}

#### Run comparisons, build + export figures ####

## Only these two comparisons have any robust-significant genera.
comparisons <- list(
  list(coral_genus = "Acropora", group_a = "Healthy_NoLesion", group_b = "Diseased_Lesion", title = "AH vs DL"),
  list(coral_genus = "Acropora", group_a = "Diseased_NoLesion", group_b = "Diseased_Lesion", title = "DU vs DL")
)

cd_list <- list()
for (spec in comparisons) {
  cmp <- run_comparison(ps_clean, ps_relab_src, spec$coral_genus, spec$group_a, spec$group_b)
  cat(sprintf("%s %s vs %s: %d tested, %d significant, %d robust\n",
              spec$coral_genus, spec$group_a, spec$group_b, cmp$n_tested, cmp$n_sig, cmp$n_robust))
  if (cmp$n_robust == 0) next
  cd_list[[spec$title]] <- build_comparison_data(cmp, spec$title)
}

## Shared x-axis (lfc) limits and bar width/height-per-row across ALL panels
all_lfc <- bind_rows(lapply(cd_list, function(cd) cd$bar_df %>% select(lfc, se)))
lfc_limits <- c(min(all_lfc$lfc - all_lfc$se), max(all_lfc$lfc + all_lfc$se)) * 1.05
bar_width <- 0.7
row_height <- 0.22    # inches per genus row, shared by bar + bubble panels alike
strip_height <- 0.9   # inches reserved for the month-shape strip + its axis text

## Render, save individually, then combine bar+bubble+strip per comparison
combined_rows <- list()
row_heights_in <- numeric(0)
for (nm in names(cd_list)) {
  cd <- cd_list[[nm]]
  p_bar <- plot_lfc_bar(cd, lfc_limits, bar_width)
  p_bubble <- plot_bubble(cd)
  p_strip <- plot_month_strip(cd)

  h_main <- max(3, row_height * cd$n_taxa + 1.5)
  h_total <- h_main + strip_height
  w_bar <- 6
  w_bubble <- max(6, 0.30 * cd$n_samples + 3)

  p_bubble_full <- (p_bubble / p_strip) + plot_layout(heights = c(h_main, strip_height))

  safe_nm <- gsub("[^A-Za-z0-9]+", "_", nm)
  ggsave(file.path(out_dir, paste0("ANCOMBC2_LFC_", safe_nm, ".pdf")), p_bar, width = w_bar, height = h_main)
  ggsave(file.path(out_dir, paste0("ANCOMBC2_Bubble_", safe_nm, ".pdf")), p_bubble_full, width = w_bubble, height = h_total)
  ggsave(file.path(out_dir, paste0("ANCOMBC2_LFC_", safe_nm, ".png")), p_bar, width = w_bar, height = h_main, dpi = 300)
  ggsave(file.path(out_dir, paste0("ANCOMBC2_Bubble_", safe_nm, ".png")), p_bubble_full, width = w_bubble, height = h_total, dpi = 300)

  ## Combined-figure
  row_design <- "AB\nCD"
  combined_rows[[nm]] <- (p_bar + p_bubble + plot_spacer() + p_strip) +
    plot_layout(design = row_design, heights = c(h_main, strip_height), widths = c(1, 1.6), guides = "collect") &
    theme(legend.position = "right")
  row_heights_in <- c(row_heights_in, h_total)

  write.csv(cd$bar_df %>% select(-TaxonLabel, -FamilyColorGroup),
            file.path(out_dir, paste0("ANCOMBC2_SigGenusRobust_", safe_nm, ".csv")), row.names = FALSE)
}

#### Combined summary table (both comparisons, all 3 groups filled in) ####

## Additional export, on top of the per-comparison ANCOMBC2_SigGenusRobust_<comparison>.csv files 
if (length(cd_list) > 0) {
  acro_meta    <- as(sample_data(ps_clean), "data.frame")
  acro_samples <- rownames(acro_meta)[acro_meta$GENUS == "Acropora"]

  ## --- ASV Number by health-status group (non-rarefied, full Acropora set) ---
  ps_acro_clean <- prune_taxa(taxa_sums(prune_samples(acro_samples, ps_clean)) > 0,
                              prune_samples(acro_samples, ps_clean))
  otu_acro <- as(otu_table(ps_acro_clean), "matrix"); if (!taxa_are_rows(ps_acro_clean)) otu_acro <- t(otu_acro)
  tx_acro  <- as(tax_table(ps_acro_clean), "matrix")[, "Taxon"]
  smp_acro <- as(sample_data(ps_acro_clean), "data.frame")
  asv_by_status <- sapply(status_abbr, function(ab) {
    cols <- rownames(smp_acro)[smp_acro$STATUS_ABBR == ab]
    tapply(rowSums(otu_acro[, cols, drop = FALSE] > 0) > 0, tx_acro, sum)
  })
  asv_by_status_df <- as.data.frame(asv_by_status) %>%
    tibble::rownames_to_column("Taxon") %>%
    dplyr::rename_with(~ paste0("n_ASVs_", .x), dplyr::all_of(status_abbr))

  ## --- Prevalence / RelAb by health-status group (RELAB_USE_RAREFY basis)
  ps_acro_relab   <- prune_samples(acro_samples, ps_relab_src)
  acro_relab_long <- shushah_relabund_long(ps_acro_relab, level = "Genus")
  acro_relab_meta <- as(sample_data(ps_acro_relab), "data.frame") %>%
    transmute(SampleID, STATUS_ABBR)
  n_in_group <- table(acro_relab_meta$STATUS_ABBR)   # Acropora sample count per status: 9 / 16 / 16

  acro_relab_wide <- shushah_group_summary(acro_relab_long, acro_relab_meta, group_cols = "STATUS_ABBR") %>%
    transmute(Taxon, STATUS_ABBR,
              Prevalence = n_present / n_samples,
              mean_RelAb = mean_RelativeAbundance, min_RelAb = min_RelativeAbundance,
              max_RelAb = max_RelativeAbundance) %>%
    tidyr::pivot_wider(names_from = STATUS_ABBR,
                       values_from = c(Prevalence, mean_RelAb, min_RelAb, max_RelAb),
                       names_glue = "{.value}_{STATUS_ABBR}")

  acro_group_stats <- dplyr::full_join(asv_by_status_df, acro_relab_wide, by = "Taxon")

  ## Taxon label: "<Genus> (<Family>)" if genus-level classified,
  ## else "Unclassified <deepest classified rank's name>"
  make_manuscript_label <- function(df) {
    ranks_deep_to_shallow <- rev(tax_ranks)
    label <- rep(NA_character_, nrow(df)); matched_rank <- rep(NA_character_, nrow(df))
    for (r in ranks_deep_to_shallow) {
      todo <- is.na(label)
      hit <- todo & !is.na(df[[r]]) & df[[r]] != "Unclassified"
      label[hit] <- df[[r]][hit]; matched_rank[hit] <- r
    }
    ifelse(matched_rank == "Genus", paste0(label, " (", df$Family, ")"), paste0("Unclassified ", label))
  }

  ## "mean (min-max)" as PERCENT, 2 decimal places, no "%" symbol
  fmt_num <- function(x) ifelse(is.na(x), NA_character_, sprintf("%.2f", x * 100))
  fmt_relab <- function(mean_v, min_v, max_v) sprintf("%s (%s-%s)", fmt_num(mean_v), fmt_num(min_v), fmt_num(max_v))

  combined_summary <- dplyr::bind_rows(lapply(names(cd_list), function(nm) {
    cd_list[[nm]]$bar_df %>%
      transmute(Taxon, Taxon_Label = make_manuscript_label(.), Group_Contrast = nm,
                logFC = lfc, SE = se, W_statistic = w_statistic, p_adj = q, n_ASVs_total = n_ASVs)
  })) %>%
    dplyr::left_join(acro_group_stats, by = "Taxon") %>%
    dplyr::mutate(
      RelAb_AH = fmt_relab(mean_RelAb_AH, min_RelAb_AH, max_RelAb_AH),
      RelAb_DU = fmt_relab(mean_RelAb_DU, min_RelAb_DU, max_RelAb_DU),
      RelAb_DL = fmt_relab(mean_RelAb_DL, min_RelAb_DL, max_RelAb_DL)
    ) %>%
    dplyr::transmute(
      Taxon = Taxon_Label, Taxon_ID = Taxon, Group_Contrast, logFC, SE, W_statistic, p_adj,
      ASV_Number_AH = n_ASVs_AH, Prevalence_AH,
      ASV_Number_DU = n_ASVs_DU, Prevalence_DU,
      ASV_Number_DL = n_ASVs_DL, Prevalence_DL,
      RelAb_AH, RelAb_DU, RelAb_DL
    )

  n_missing <- sum(!stats::complete.cases(
    combined_summary %>% dplyr::select(dplyr::starts_with("ASV_Number"), dplyr::starts_with("Prevalence"))
  ))
  if (n_missing > 0) {
    warning(n_missing, " row(s) in the combined summary table are missing a group's",
            " ASV/prevalence stats - check acro_group_stats coverage.", call. = FALSE)
  }

  ## Column headers: sample-group size spelled out in the Prevalence header
  final_names <- c(
    Taxon = "Taxon", Taxon_ID = "Taxon_ID", Group_Contrast = "Group Contrast",
    logFC = "logFC", SE = "SE", W_statistic = "W_statistic", p_adj = "p-adj",
    ASV_Number_AH = "ASV Number AH",
    Prevalence_AH = sprintf("Per-group Prevalence AH (out of %d samples)", n_in_group[["AH"]]),
    ASV_Number_DU = "ASV Number DU",
    Prevalence_DU = sprintf("Per-group Prevalence DU (out of %d samples)", n_in_group[["DU"]]),
    ASV_Number_DL = "ASV Number DL",
    Prevalence_DL = sprintf("Per-group Prevalence DL (out of %d samples)", n_in_group[["DL"]]),
    RelAb_AH = "mean RelAb AH (min-max)",
    RelAb_DU = "mean RelAb DU (min-max)",
    RelAb_DL = "mean RelAb DL (min-max)"
  )
  names(combined_summary) <- unname(final_names[names(combined_summary)])

  write.csv(combined_summary, file.path(out_dir, "ANCOMBC2_SigGenusRobust_Combined.csv"), row.names = FALSE)
  writeLines(c(
    "ANCOMBC2_SigGenusRobust_Combined.csv - column notes",
    "====================================================",
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "Source script: ShushahTotal_Microbiome_ANCOMBC2_forPublishing.R"
  ), file.path(out_dir, "ANCOMBC2_SigGenusRobust_Combined_METHODS.txt"))
  cat("\nCombined manuscript summary table: ", nrow(combined_summary), " rows -> ",
      "ANCOMBC2_SigGenusRobust_Combined.csv\n", sep = "")
}

## Full combined figure: one row per comparison, bar + bubble(+strip) columns
if (length(combined_rows) > 0) {
  full_figure <- wrap_plots(combined_rows, ncol = 1) +
    plot_layout(heights = row_heights_in) +
    plot_annotation(tag_levels = "a")
  fig_h <- sum(row_heights_in)
  fig_w <- max(sapply(cd_list, function(cd) max(6, 0.30 * cd$n_samples + 3))) + 6
  ggsave(file.path(out_dir, "ANCOMBC2_Combined_Figure.pdf"), full_figure, width = fig_w, height = fig_h, limitsize = FALSE)
  ggsave(file.path(out_dir, "ANCOMBC2_Combined_Figure.png"), full_figure, width = fig_w, height = fig_h, dpi = 300, limitsize = FALSE)
}

cat("\nDone. Outputs in", out_dir, "\n")
