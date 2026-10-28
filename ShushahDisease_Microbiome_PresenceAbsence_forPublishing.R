#### Presence/Absence Analysis - Shushah Disease Microbiome ####
#### Libraries ####

library(tidyverse)   # data wrangling / ggplot2
library(phyloseq)    # microbiome data object + manipulation
library(microbiome)  # convenience helpers
library(openxlsx)    # read the approved family-color palette

# Shared relative-abundance helpers to compute relative abundances consistently across codes
source("./shushah_relabund_helpers.R") #Edit with correct path

#### Configuration ####

## Working directory (contains the raw ASV/taxonomy/metadata files)
setwd(
  "./ShushahDisease_Microbiome_PresenceAbsence" #Edit with correct path
)
dir.create("Results", showWarnings = FALSE)

## Taxonomic resolution for the whole analysis
analysis_level <- "Genus"

## Full taxonomic rank order in the taxonomy file (Kingdom -> Genus), used to
## build disambiguated taxon labels when agglomerating (see glom_to_level()).
tax_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

## Rarefaction depth + seed. Fixed at 42800 to match the depth used
## throughout the rest of the project's ANCOM-BC2/composition analyses
rarefy_depth <- 42800
rarefy_seed <- 69

RELAB_USE_RAREFY <- TRUE
## TRUE  → Mean/Min/Max relative abundance values and plots are generated on the rarefied datasets
## FALSE → Uses cleaned, non rarefied dataset for the values above
## Note: presence/absence and exclusivity are always computed on the rarefied object, regardless of this switch

#### Import raw data ####

taxfile_PA <- read.csv(
  "ShushahMicrobiome_Total_ASVs_taxonomy-no-contam.tsv",
  header = TRUE, row.names = 1, sep = "\t"
)

mapfile_PA <- read.csv(
  "ShushahMicrobiome_Total_Metadata.csv",
  header = TRUE, row.names = 1
)

ASVs_PA <- read.csv(
  "ShushahMicrobiome_Total_ASVs_counts-no-contam.tsv",
  header = TRUE, row.names = 1, sep = "\t"
)

## Matrix conversions expected by phyloseq
taxfile_PA <- as.matrix(taxfile_PA)
ASVs_PA <- as.matrix(ASVs_PA)

## Tidy metadata: add SampleID + a single 3-level health-status factor
mapfile_PA <- data.frame(mapfile_PA)
mapfile_PA$SampleID <- rownames(mapfile_PA)

mapfile_PA <- mapfile_PA %>%
  unite(col = "STATUS", COLONY_STATUS, FRAGMENT_STATUS, sep = "_", remove = FALSE)

## Reorder and add short labels
status_levels <- c("Healthy_NoLesion", "Diseased_NoLesion", "Diseased_Lesion")
status_abbr   <- c("AH", "DU", "DL")

mapfile_PA$STATUS <- factor(mapfile_PA$STATUS, levels = status_levels)
mapfile_PA$STATUS_ABBR <- factor(
  dplyr::recode(mapfile_PA$STATUS,
         "Healthy_NoLesion" = "AH",
         "Diseased_NoLesion" = "DU",
         "Diseased_Lesion" = "DL"),
  levels = status_abbr
)
mapfile_PA$MONTH <- factor(mapfile_PA$MONTH, levels = c("March", "July"))

#### Build phyloseq object ####
tree_PA <- phytools::midpoint.root(
  phyloseq::read_tree("ShushahMicrobiome_Total_ASVs-no-contam_tree.nwk")
)

PA_phyloseq <- phyloseq(
  otu_table(ASVs_PA, taxa_are_rows = TRUE),
  tax_table(taxfile_PA),
  sample_data(mapfile_PA),
  phy_tree(tree_PA)
)
PA_phyloseq

## Organelle/non-bacterial + global-singleton filter. 
## Shared helper, identical across all three project scripts
PA_phyloseq_clean <- shushah_clean_phyloseq(PA_phyloseq)
PA_phyloseq_clean

## Number of ASVs that originally belong to each disambiguated genus ("Taxon"),
## counted across the FULL cleaned (non-rarefied, all 65 samples) dataset
if (!identical(analysis_level, "ASV")) {
  taxon_asv_n_PA <- as.data.frame(
    table(Taxon = as(tax_table(shushah_add_taxon_column(PA_phyloseq_clean,
                                                        level = analysis_level)),
                      "matrix")[, "Taxon"]),
    responseName = "n_ASVs", stringsAsFactors = FALSE
  )
  taxon_asv_n_PA$Taxon  <- as.character(taxon_asv_n_PA$Taxon)
  taxon_asv_n_PA$n_ASVs <- as.integer(taxon_asv_n_PA$n_ASVs)
} else {
  taxon_asv_n_PA <- data.frame(Taxon = character(0), n_ASVs = integer(0))
}

#### Rarefaction ####

cat("Rarefying all samples to", rarefy_depth, "reads.\n")

PA_phyloseq_rar <- shushah_rarefy(PA_phyloseq_clean, depth = rarefy_depth, seed = rarefy_seed)
PA_phyloseq_rar

## Relative-abundance source object (see RELAB_USE_RAREFY in Configuration).
PA_relabund_src <- if (isTRUE(RELAB_USE_RAREFY)) PA_phyloseq_rar else PA_phyloseq_clean

#### Subset by coral genus ####

PA_Acropora <- PA_phyloseq_rar %>%
  subset_samples(GENUS == "Acropora") %>%
  { prune_taxa(taxa_sums(.) > 0, .) }

PA_Pocillopora <- PA_phyloseq_rar %>%
  subset_samples(GENUS == "Pocillopora") %>%
  { prune_taxa(taxa_sums(.) > 0, .) }

PA_Acropora
PA_Pocillopora

#### Agglomerate to the configured analysis_level ####

## Assign a disambiguated taxon name first to avoid collapsing "Unclassified" genera into a single entity
glom_to_level <- function(ps, level) {
  if (identical(level, "ASV")) return(ps)

  rank_idx <- which(tax_ranks == level)
  if (length(rank_idx) == 0) stop("analysis_level must be \"ASV\" or one of: ", paste(tax_ranks, collapse = ", "))
  upstream_ranks <- tax_ranks[seq_len(rank_idx - 1)]

  tt_df <- as.data.frame(as(tax_table(ps), "matrix"), stringsAsFactors = FALSE)
  lineage <- if (length(upstream_ranks) > 0) {
    do.call(paste, c(tt_df[upstream_ranks], sep = "_"))
  } else {
    rep("Root", nrow(tt_df))
  }
  tt_df$Taxon <- ifelse(
    !is.na(tt_df[[level]]) & tt_df[[level]] != "Unclassified",
    tt_df[[level]],
    paste(lineage, "Unclassified", sep = "_")
  )

  new_tax <- as.matrix(tt_df)
  rownames(new_tax) <- taxa_names(ps)
  tax_table(ps) <- tax_table(new_tax)
  tax_glom(ps, taxrank = "Taxon", NArm = FALSE)
}

cat("Agglomerating to analysis_level =", analysis_level, "\n")
PA_Acropora <- glom_to_level(PA_Acropora, analysis_level)
PA_Pocillopora <- glom_to_level(PA_Pocillopora, analysis_level)

## Parallel genus subsets built from the relative-abundance source object
## (rarefied or non-rarefied, per RELAB_USE_RAREFY), agglomerated identically.
## When RELAB_USE_RAREFY = TRUE these are identical to PA_Acropora / PA_Pocillopora.
if (isTRUE(RELAB_USE_RAREFY)) {
  PA_Acropora_relab    <- PA_Acropora
  PA_Pocillopora_relab <- PA_Pocillopora
} else {
  PA_Acropora_relab <- PA_relabund_src %>%
    subset_samples(GENUS == "Acropora") %>%
    { prune_taxa(taxa_sums(.) > 0, .) } %>%
    glom_to_level(analysis_level)
  PA_Pocillopora_relab <- PA_relabund_src %>%
    subset_samples(GENUS == "Pocillopora") %>%
    { prune_taxa(taxa_sums(.) > 0, .) } %>%
    glom_to_level(analysis_level)
}

PA_Acropora
PA_Pocillopora

#### Helper functions ####

tax_df <- function(ps) {
  as.data.frame(as(tax_table(ps), "matrix")) %>%
    rownames_to_column("TaxonID")
}

presence_matrix <- function(ps) {
  as(otu_table(ps), "matrix") > 0
}

relabund_matrix <- function(ps) {
  m <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) m <- t(m)
  m <- sweep(m, 2, colSums(m), "/")
  tt <- as(tax_table(ps), "matrix")
  if ("Taxon" %in% colnames(tt)) rownames(m) <- tt[, "Taxon"]
  m
}

## --- Group-exclusive taxa ----------------------------------------------
## A taxon is "group-exclusive" if it is present in samples of only one STATUS group
## Always uses the RAREFIED object.
## `ps_relab` supplies relative-abundance magnitudes, computed over EVERY sample of
## the exclusive group, including samples where the taxon was not detected 

get_group_exclusive <- function(ps, ps_relab = ps) {
  pres  <- presence_matrix(ps)
  relab <- relabund_matrix(ps_relab)
  tx_map <- as(tax_table(ps), "matrix")[, "Taxon"]   # named by TaxonID
  smp <- sample_data(ps) %>% as("data.frame")

  groups <- levels(droplevels(smp$STATUS_ABBR))
  group_samples <- split(smp$SampleID, smp$STATUS_ABBR)[groups]
  group_n <- lengths(group_samples)

  ## taxon x group presence (present in >=1 sample of that group)
  group_pres <- sapply(groups, function(g) {
    cols <- group_samples[[g]]
    if (length(cols) == 1) pres[, cols] else rowSums(pres[, cols, drop = FALSE]) > 0
  })
  rownames(group_pres) <- rownames(pres)

  n_groups_present <- rowSums(group_pres)
  exclusive_ids <- rownames(group_pres)[n_groups_present == 1]

  if (length(exclusive_ids) == 0) {
    return(tibble(
      TaxonID = character(), ExclusiveGroup = character(),
      N_SamplesPresent = integer(), N_SamplesInGroup = integer(), Prevalence = double(),
      MeanRelativeAbundance = double(), MinRelativeAbundance = double(),
      MaxRelativeAbundance = double(), n_ASVs = integer(), SampleIDs = character()
    ))
  }

  map_dfr(exclusive_ids, function(id) {
    g <- groups[group_pres[id, ]]
    cols <- group_samples[[g]]                      # ALL samples of the exclusive group
    present_cols <- cols[pres[id, cols]]            # detected subset (for prevalence only)
    tx <- tx_map[[id]]
    ra <- if (!is.na(tx) && tx %in% rownames(relab)) relab[tx, cols] else rep(NA_real_, length(cols))
    tibble(
      TaxonID = id,
      ExclusiveGroup = g,
      N_SamplesPresent = length(present_cols),
      N_SamplesInGroup = group_n[[g]],
      Prevalence = length(present_cols) / group_n[[g]],
      MeanRelativeAbundance = mean(ra),
      MinRelativeAbundance  = min(ra),
      MaxRelativeAbundance  = max(ra),
      SampleIDs = paste(present_cols, collapse = "; ")
    )
  }) %>%
    left_join(tax_df(ps), by = "TaxonID") %>%
    { if ("Taxon" %in% names(.)) left_join(., taxon_asv_n_PA, by = "Taxon") else . } %>%
    relocate(any_of(c("TaxonID", "ExclusiveGroup", "N_SamplesPresent", "N_SamplesInGroup",
                      "Prevalence", "MeanRelativeAbundance", "MinRelativeAbundance",
                      "MaxRelativeAbundance", "n_ASVs", "SampleIDs")))
}

#### Group-exclusive taxa, per genus ####

PA_groupExcl_Acropora <- get_group_exclusive(PA_Acropora, PA_Acropora_relab)
PA_groupExcl_Pocillopora <- get_group_exclusive(PA_Pocillopora, PA_Pocillopora_relab)

cat("\nAcropora: ", nrow(PA_groupExcl_Acropora), " group-exclusive taxa (unfiltered)\n", sep = "")
cat("Pocillopora: ", nrow(PA_groupExcl_Pocillopora), " group-exclusive taxa (unfiltered)\n", sep = "")

#### Export CSVs ####

write.csv(PA_groupExcl_Acropora, "Results/PresenceAbsence_Acropora_Combined.csv", row.names = FALSE)
write.csv(PA_groupExcl_Pocillopora, "Results/PresenceAbsence_Pocillopora_Combined.csv", row.names = FALSE)

cat("\nCSV tables written:\n",
    " - Results/PresenceAbsence_Acropora_Combined.csv\n",
    " - Results/PresenceAbsence_Pocillopora_Combined.csv\n")

#### Prevalence filter ####
## A taxon is dropped unless it is present in at least this fraction of its
## own exclusive group's samples (e.g. 0.80 -> at least 8/10 samples).
## Set to Null to disable filtering
## P.S. filtering is applied after the full table above, before plotting/export

prevalence_threshold <- 0.8

apply_prevalence_filter <- function(group_excl, threshold) {
  if (is.null(threshold)) return(group_excl)
  group_excl %>% filter(Prevalence >= threshold)
}

PA_groupExcl_Acropora_filtered <- apply_prevalence_filter(PA_groupExcl_Acropora, prevalence_threshold)
PA_groupExcl_Pocillopora_filtered <- apply_prevalence_filter(PA_groupExcl_Pocillopora, prevalence_threshold)

cat("\nAfter prevalence filtering (>=", prevalence_threshold, "):\n", sep = "")
cat("Acropora: ", nrow(PA_groupExcl_Acropora_filtered), " of ", nrow(PA_groupExcl_Acropora), " taxa retained\n", sep = "")
cat("Pocillopora: ", nrow(PA_groupExcl_Pocillopora_filtered), " of ", nrow(PA_groupExcl_Pocillopora), " taxa retained\n", sep = "")

## Sub-folder named after the threshold in use, e.g. "Results_prevalence_0.80"
prevalence_dir <- file.path("Results", sprintf(
  "Results_prevalence_%s",
  if (is.null(prevalence_threshold)) "none" else sprintf("%.2f", prevalence_threshold)
))
dir.create(prevalence_dir, showWarnings = FALSE)

write.csv(PA_groupExcl_Acropora_filtered, file.path(prevalence_dir, "PresenceAbsence_Acropora_Combined.csv"), row.names = FALSE)
write.csv(PA_groupExcl_Pocillopora_filtered, file.path(prevalence_dir, "PresenceAbsence_Pocillopora_Combined.csv"), row.names = FALSE)

cat("\nFiltered CSV tables written to ", prevalence_dir, "/\n", sep = "")

#### Relative abundance exports ####
## Independent of the group-exclusivity analysis and its prevalence filter

relab_export_dir <- file.path("Results", "relative_abundance")
dir.create(relab_export_dir, showWarnings = FALSE)

relab_sample_meta <- sample_data(PA_relabund_src) %>% as("data.frame") %>%
  transmute(SampleID, GENUS = as.character(GENUS), STATUS_ABBR)

for (lvl in c("Family", "Genus")) {
  relab_long  <- shushah_relabund_long(PA_relabund_src, level = lvl)
  relab_summ  <- shushah_group_summary(relab_long, relab_sample_meta,
                                       group_cols = c("GENUS", "STATUS_ABBR"))
  relab_top19 <- shushah_topN_union_summary(relab_long, relab_sample_meta,
                                            group_cols = c("GENUS", "STATUS_ABBR"),
                                            subset_col = "GENUS", top_n = 19)
  write.csv(relab_long,  file.path(relab_export_dir, paste0(lvl, "_relabundance_persample.csv")),
            row.names = FALSE)
  write.csv(relab_summ,  file.path(relab_export_dir, paste0(lvl, "_relabundance_group_summary.csv")),
            row.names = FALSE)
  write.csv(relab_top19, file.path(relab_export_dir, paste0(lvl, "_relabundance_top19_group_summary.csv")),
            row.names = FALSE)
  cat("  ", lvl, ": ", nrow(relab_long), " per-sample rows, ", nrow(relab_summ),
      " group rows, ", length(attr(relab_top19, "named_taxa")),
      " named + Others (top-19 union)\n", sep = "")
}

shushah_write_relab_methods(
  file.path(relab_export_dir, "RELAB_METHODS.txt"),
  use_rarefied = isTRUE(RELAB_USE_RAREFY),
  group_desc   = paste("coral genus x health-status (AH = Healthy/NoLesion,",
                       "DU = Diseased/NoLesion, DL = Diseased/Lesion); March + July pooled"),
  extra_lines  = c(
    "Source script: ShushahDisease_Microbiome_PresenceAbsence_forPublishing.R",
    "These tables are relative-abundance summaries only - NOT the group-exclusive",
    "taxa analysis. Presence/absence detection and exclusivity always use the",
    "rarefied object regardless of RELAB_USE_RAREFY.",
    "In the PresenceAbsence_<genus>_Combined.csv exclusive-taxa tables, the",
    "MeanRelativeAbundance / MinRelativeAbundance / MaxRelativeAbundance columns",
    "are taken over EVERY sample of the exclusive group (non-detections = 0),",
    "not only the samples where the taxon was detected."
  )
)
cat("Relative abundance tables written to ", relab_export_dir, "/\n", sep = "")

#### Bubble plots ####

## size = relative abundance, on a fixed scale shared by both genus plots
## color = bacterial Family; families NOT in the palette below collapse 
##         into its own "Others" entry 

## Readable per-row label built from the deepest classified rank available
make_taxon_label <- function(df) {
  ranks_deep_to_shallow <- rev(tax_ranks)  # Genus, Family, Order, Class, Phylum, Kingdom
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

## Family-color palette shared across codes with the project's family-level relative-abundance figure
family_palette_file <- "./Color codes families SVG.xlsx" ## Edit with correct path
approved_colors <- read.xlsx(family_palette_file, sheet = 1) %>%
  mutate(HexColor = paste0("#", substr(Color, 1, 6))) %>%
  { setNames(.$HexColor, .$Family) }

bubble_plot_data <- function(ps, group_excl, ps_relab = ps) {
  if (nrow(group_excl) == 0) return(NULL)

  relab <- relabund_matrix(ps_relab)   # rows keyed by disambiguated "Taxon"
  smp <- sample_data(ps) %>% as("data.frame") %>%
    select(SampleID, STATUS_ABBR, MONTH) %>%
    mutate(STATUS_ABBR = factor(STATUS_ABBR, levels = status_abbr), MONTH = factor(MONTH, levels = c("March", "July")))

  surviving_groups <- unique(group_excl$ExclusiveGroup)
  smp <- smp %>% filter(STATUS_ABBR %in% surviving_groups) %>% mutate(STATUS_ABBR = droplevels(STATUS_ABBR))

  group_excl <- group_excl %>%
    mutate(
      FamilyColorGroup = ifelse(Family %in% names(approved_colors), Family, "Others"),
      BaseLabel = make_taxon_label(.),
      TaxonLabel = ifelse(
        FamilyColorGroup == "Others" & !str_detect(BaseLabel, fixed(Family)),
        paste0(BaseLabel, " (", Family, ")"),
        BaseLabel
      )
    ) %>%
    mutate(TaxonLabel = make.unique(TaxonLabel, sep = " #"))

  row_order_df <- group_excl %>%
    arrange(desc(MeanRelativeAbundance))
  taxon_order <- rev(row_order_df$TaxonLabel)

  sample_order <- smp %>%
    arrange(MONTH, STATUS_ABBR, SampleID) %>%
    mutate(SampleLabel = paste(SampleID, STATUS_ABBR, MONTH, sep = "\n"))

  bubble_data <- map_dfr(seq_len(nrow(group_excl)), function(i) {
    ids <- intersect(str_split(group_excl$SampleIDs[i], "; ")[[1]], sample_order$SampleID)
    if (length(ids) == 0) return(NULL)
    tibble(TaxonID = group_excl$TaxonID[i], Taxon = group_excl$Taxon[i], SampleID = ids)
  }) %>%
    mutate(RelativeAbundance = relab[cbind(Taxon, SampleID)]) %>%
    left_join(group_excl %>% select(TaxonID, TaxonLabel, FamilyColorGroup), by = "TaxonID") %>%
    left_join(sample_order %>% select(SampleID, SampleLabel), by = "SampleID") %>%
    mutate(TaxonLabel = factor(TaxonLabel, levels = taxon_order), SampleLabel = factor(SampleLabel, levels = sample_order$SampleLabel))

  color_group_levels <- row_order_df %>%
    mutate(ColorSortKey = case_when(
      FamilyColorGroup == "Unclassified" ~ "zzz2_Unclassified",
      FamilyColorGroup == "Others" ~ "zzz1_Others",
      TRUE ~ FamilyColorGroup
    )) %>%
    distinct(FamilyColorGroup, ColorSortKey) %>%
    arrange(ColorSortKey) %>%
    pull(FamilyColorGroup)

  list(data = bubble_data, color_group_levels = color_group_levels,
       family_palette = approved_colors[color_group_levels],
       n_taxa = length(taxon_order), n_samples = nrow(sample_order))
}

## Fixed bubble-size scale shared across BOTH genus plots.
make_size_scale <- function(global_max) {
  breaks <- scales::extended_breaks(n = 5)(c(0, global_max))
  breaks <- breaks[breaks > 0 & breaks <= global_max]
  scale_size_area(
    max_size = 8, limits = c(0, global_max), breaks = breaks,
    labels = scales::percent_format(accuracy = 0.1), name = "Relative\nabundance"
  )
}

plot_bubble <- function(bd, genus_name, size_scale) {
  bd$data <- bd$data %>% mutate(FamilyColorGroup = factor(FamilyColorGroup, levels = bd$color_group_levels))
  ggplot(bd$data, aes(x = SampleLabel, y = TaxonLabel, size = RelativeAbundance, fill = FamilyColorGroup)) +
    geom_point(shape = 21, color = "grey30", stroke = 0.15, alpha = 0.9) +
    size_scale +
    scale_fill_manual(values = bd$family_palette, name = "Family", drop = FALSE) +
    ## Samples with zero bubbles must still show as an empty column
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    guides(fill = guide_legend(override.aes = list(size = 4), ncol = 1)) +
    labs(
      x = NULL, y = NULL,
      title = str_wrap(paste0(
        genus_name, ": genera exclusive to one health-status group (Combined months, prevalence >= ",
        scales::percent(prevalence_threshold), ")"
      ), 55)
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(size = 6, lineheight = 0.8),
      axis.text.y = element_text(size = 8),
      panel.grid.major = element_line(color = "grey93"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 11)
    )
}

bd_list <- list()
for (genus_name in c("Acropora", "Pocillopora")) {
  ps       <- if (genus_name == "Acropora") PA_Acropora else PA_Pocillopora
  ps_relab <- if (genus_name == "Acropora") PA_Acropora_relab else PA_Pocillopora_relab
  group_excl <- if (genus_name == "Acropora") PA_groupExcl_Acropora_filtered else PA_groupExcl_Pocillopora_filtered

  bd <- bubble_plot_data(ps, group_excl, ps_relab)
  if (is.null(bd)) {
    cat("\nNo group-exclusive taxa for", genus_name, "after prevalence filtering - skipping bubble plot.\n")
    next
  }
  bd_list[[genus_name]] <- bd
}

global_max_relab <- max(unlist(lapply(bd_list, function(bd) bd$data$RelativeAbundance)), na.rm = TRUE)
size_scale <- make_size_scale(global_max_relab)
cat("\nShared bubble-size scale: 0 to", scales::percent(global_max_relab, accuracy = 0.1), "\n")

for (genus_name in names(bd_list)) {
  bd <- bd_list[[genus_name]]
  w <- max(7, 0.30 * bd$n_samples + 4)
  h <- max(5, 0.20 * bd$n_taxa + 2.5)

  p <- plot_bubble(bd, genus_name, size_scale)

  out_base <- file.path(prevalence_dir, sprintf("PresenceAbsence_Bubble_%s", genus_name))
  ggsave(paste0(out_base, ".pdf"), p, width = w, height = h, limitsize = FALSE)
  ggsave(paste0(out_base, ".png"), p, width = w, height = h, dpi = 200, limitsize = FALSE)

  cat("\n", genus_name, " bubble plot: ", bd$n_taxa, " taxa x ", bd$n_samples,
      " samples, ", length(bd$color_group_levels), " color groups -> ", out_base, ".pdf/.png\n", sep = "")
}
