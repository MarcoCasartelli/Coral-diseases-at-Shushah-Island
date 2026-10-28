#### Coral diseaseas in Shushah Island - Microbiome composition pipeline ####

# Edit SECTION 1 to SECTION 8 for the desired analysis
# Each section and subsection can be enabled/disabled with its top-level flag (TRUE/FALSE)

#### Libraries ####

required_packages <- c(
  "tidyverse",
  "phyloseq",
  "phytools",
  "vegan",
  "dunn.test",
  "paletteer",
  "car",      
  "ggpubr"    
)

invisible(lapply(required_packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package not installed: ", pkg, call. = FALSE)
  }
  library(pkg, character.only = TRUE)
}))

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

# Shared relative-abundance helpers to compute relative abundances consistently across codes
source("/Users/casartm/Desktop/KAUST PhD/Projects/Coral Disease - Shushah/Microbiome Analysis/ShushahDisease_RelAb_Helpers/shushah_relabund_helpers.R")

# ==============================================================================
# SECTION 1 — PATHS & PROJECT
# ==============================================================================

PROJECT_NAME <- "ShushahTotal_Disease"

# Folder containing ASV taxonomy, ASV count, metadata, and (optionally) tree files. 
# Change to "." if those files are in the project root.
DATA_DIR <- "/Users/casartm/Desktop/KAUST PhD/Projects/Coral Disease - Shushah/Microbiome Analysis/ShushahDisease_MicrobiomeComposition"

OUTPUT_DIR <- "/Users/casartm/Desktop/KAUST PhD/Projects/Coral Disease - Shushah/Microbiome Analysis/ShushahDisease_MicrobiomeComposition"

# File names inside DATA_DIR
FILE_TAXONOMY <- "ShushahMicrobiome_Total_ASVs_taxonomy-no-contam.tsv"
FILE_METADATA <- "ShushahMicrobiome_Total_Metadata.csv"
FILE_COUNTS   <- "ShushahMicrobiome_Total_ASVs_counts-no-contam.tsv"

# Tree is OPTIONAL — set to NULL (or "") if you don't have a tree file.
# When no tree is available, the phyloseq object is built without one, and the UniFrac 
# distance methods in SECTION 6c are automatically skipped at run time with a console message.
FILE_TREE <- "ShushahMicrobiome_Total_ASVs-no-contam_tree.nwk"

# Column in metadata that uniquely identifies each sample
SAMPLE_ID_COLUMN <- "Sample_ID"

# ==============================================================================
# SECTION 2 — METADATA & GROUPING
# ==============================================================================

# The grouping variable used throughout the pipeline.
GROUP_VAR <- "GROUP"

# Columns to paste together to create GROUP_VAR (set to c() to skip)
GROUP_FROM <- c("GENUS", "COLONY_STATUS", "FRAGMENT_STATUS")
GROUP_SEP  <- "_"

# Ordered factor levels for GROUP_VAR (set to NULL for unordered)
GROUP_LEVELS <- c(
  "Acropora_Healthy_NoLesion",
  "Acropora_Diseased_NoLesion",
  "Acropora_Diseased_Lesion",
  "Pocillopora_Healthy_NoLesion",
  "Pocillopora_Diseased_NoLesion",
  "Pocillopora_Diseased_Lesion"
)

# Ordered factor levels for MONTH (set to NULL for unordered)
SAMPLE_ORIGIN_LEVELS <- c(
  "March",
  "July"
)

# Which ordered levels apply to which metadata column
VARIABLE_LEVELS <- list(
  "GROUP"         = GROUP_LEVELS,
  "MONTH" = SAMPLE_ORIGIN_LEVELS
)

# Combined column pasting SAMPLE_ID_COLUMN and GROUP_VAR together 
# (e.g. TAXONOMY_SAMPLE_VAR in SECTION 7). 
# Set to NULL to skip creating it.
SAMPLE_GROUP_LABEL_COLUMN <- "Sample_ID_Group"
SAMPLE_GROUP_LABEL_SEP    <- "_"

# ==============================================================================
# SECTION 3 — RAREFACTION
# ==============================================================================

RAREFACTION_ENABLED     <- TRUE
RAREFACTION_SAMPLE_SIZE <- 42800
RAREFACTION_SEED        <- 69

# --- 3a. Sequencing depth diagnostics (optional) ---
# Before rarefying, export a CSV of read counts per sample and a rarefaction curve 
# (ASVs discovered vs. sequencing depth, one line per sample)
RAREFACTION_DIAGNOSTICS <- FALSE
RAREFACTION_CURVE_STEP  <- 10000   # step size (in reads) used by rarecurve()

# ==============================================================================
# SECTION 4 — SAMPLE SUBSETS
# ==============================================================================

# Add or edit subsets here. Any metadata column can be used in include/exclude.
#   include: keep rows where column value is in values list
#   exclude: drop rows where column value is in values list

SUBSETS <- list(
  Overall = list(
    title   = "All Samples",
    include = list(),
    exclude = list()
  # ),
  # Acropora = list(
  #   title   = "Acropora",
  #   include = list(GENUS = "Acropora"),
  #   exclude = list()
  # ),
  # Pocillopora = list(
  #   title   = "Pocillopora",
  #   include = list(GENUS = "Pocillopora"),
  #   exclude = list()
  # ),
  # March = list(
  #   title   = "March",
  #   include = list(MONTH = "March"),
  #   exclude = list()
  # ),
  # July = list(
  #   title   = "July",
  #   include = list(MONTH = "July"),
  #   exclude = list()
  # ),
  # Acropora_March = list(
  #   title   = "Acropora - March",
  #   include = list(GENUS = "Acropora"),
  #   exclude = list(MONTH = "July")
  # ),
  # Acropora_July = list(
  #   title   = "Acropora - July",
  #   include = list(GENUS = "Acropora"),
  #   exclude = list(MONTH = "March")
  # ),
  # Pocillopora_March = list(
  #   title   = "Pocillopora - March",
  #   include = list(GENUS = "Pocillopora"),
  #   exclude = list(MONTH = "July")
  # ),
  # Pocillopora_July = list(
  #   title   = "Pocillopora - July",
  #   include = list(GENUS = "Pocillopora"),
  #   exclude = list(MONTH = "March")
  # ),
  # Acropora_Species = list(
  #   title   = "Acropora_Species",
  #   include = list(GENUS = "Acropora"),
  #   exclude = list(SPECIES = c("A. cf. cytherea", "A. maryae"),
  #                  FRAGMENT_STATUS = "Lesion")
  # ),
  # Acropora_Species_March = list(
  #   title   = "Acropora_Species",
  #   include = list(GENUS = "Acropora",
  #                  MONTH = "March"),
  #   exclude = list(SPECIES = c("A. cf. cytherea", "A. maryae"),
  #                  FRAGMENT_STATUS = "Lesion")
  # ),
  # Diseased_Colonies = list(
  #   title   = "Diseased Colonies (Lesion vs NoLesion, paired by colony)",
  #   include = list(COLONY_STATUS = "Diseased"),
  #   exclude = list()
  # ),
  # Diseased_Acropora = list(
  #   title   = "Diseased Colonies (Lesion vs NoLesion, paired by colony)",
  #   include = list(GENUS = "Acropora",
  #                  COLONY_STATUS = "Diseased"),
  #   exclude = list()
  # ),
  # Diseased_Pocillopora = list(
  #   title   = "Diseased Colonies (Lesion vs NoLesion, paired by colony)",
  #   include = list(GENUS = "Pocillopora",
  #                  COLONY_STATUS = "Diseased"),
  #   exclude = list()
  )
)

# ==============================================================================
# SECTION 5 — ALPHA DIVERSITY
# ==============================================================================

ALPHA_ENABLED      <- TRUE   # Master switch for the entire alpha diversity section

ALPHA_USE_RAREFIED <- TRUE
# TRUE  → run on the rarefied phyloseq object
# FALSE → run on the cleaned, non-rarefied phyloseq object

# --- 5a. Metrics ---
ALPHA_METRICS      <- c("Observed", "Chao1", "Shannon", "Simpson", "Evenness")
ALPHA_PLOT_METRICS <- c("Observed", "Shannon", "Simpson", "Evenness")

# Optional fixed y-axis limits per metric, shared across every subset
# NULL (default) = auto-scale each plot to its own data
# To enable, use a named list with one c(min, max) per metric you want fixed
# e.g.: ALPHA_Y_LIMITS <- list(Observed = c(0, 2000), Shannon = c(1, 7), Evenness = c(0,1))
# Metrics left out of the list still auto-scale.
ALPHA_Y_LIMITS <-NULL

# Metadata variable used for grouping in alpha tests and plots
ALPHA_GROUP_VAR <- "GROUP"

# --- 5b. Statistical testing ---
# Significance threshold used for:
#   (1) assumption tests  — groups failing p < threshold trigger non-parametric path
#   (2) omnibus test      — only run post-hoc if omnibus p < threshold
#   (3) pairwise display  — only brackets with adj. p < threshold are drawn on plots
ALPHA_SIG_THRESHOLD <- 0.05

# Multiple-comparison correction method for pairwise tests
# Dunn (non-parametric): any method accepted by p.adjust() e.g. "BH", "bonferroni"
# Tukey (parametric):    always uses its own family-wise correction internally
ALPHA_PAIRWISE_ADJUST <- "BH"

# --- 5c. Sub-feature toggles ---
ALPHA_PLOT       <- TRUE   # Generate per-metric boxplots
ALPHA_STATS      <- TRUE   # Run assumption tests + omnibus + post-hoc
ALPHA_ANNOTATE   <- TRUE   # Add significance brackets to boxplots
ALPHA_SHOW_POINTS  <- TRUE   # Overlay individual sample dots on boxplots
ALPHA_ALL_PVALUES  <- FALSE   # Show p-values for ALL pairwise comparisons, not only significant ones

# --- 5d. Optional second grouping variable (facet / point shape) ---
# When used as point shape, it reuses the shape mapping from SHAPE_BY_VARIABLE (SECTION 8b)
ALPHA_VARIABLE2_ENABLED  <- FALSE
ALPHA_VARIABLE2_COLUMN   <- "MONTH"
ALPHA_VARIABLE2_AS_FACET <- TRUE    # Facet boxplots by Variable_2
ALPHA_VARIABLE2_AS_SHAPE <- TRUE    # Use Variable_2 as point shape for the overlaid sample dots (requires ALPHA_SHOW_POINTS)

# NOTE: significance brackets (ALPHA_ANNOTATE) are computed on the whole
# dataset for test_group_var only, not per facet — if faceting is enabled,
# the same brackets are repeated identically in every panel rather than
# reflecting within-facet significance.

# ==============================================================================
# SECTION 6 — BETA DIVERSITY
# ==============================================================================

BETA_ENABLED      <- TRUE   # Master switch for the entire beta diversity section

BETA_USE_RAREFIED <- TRUE
# TRUE  → distance methods with source = "rel_abundance" draw from rarefied
#         relative abundances; source = "clean_counts" is unaffected (always
#         uses clean non-rarefied counts for Aitchison / CLR-based distances).
# FALSE → source = "rel_abundance" draws from relative abundances computed on
#         the cleaned, non-rarefied counts instead.

# --- 6a. Ordination & grouping ---
BETA_ORDINATION_METHOD <- "PCoA"
BETA_GROUP_VAR         <- "GROUP"
BETA_SAMPLE_LABEL_VAR  <- "Sample_ID"
BETA_PERMUTATIONS      <- 999
BETA_SEED              <- 69
BETA_PAIRWISE_P_ADJUST <- "BH"

# PCoA plot appearance
BETA_SHOW_SEGMENTS <- TRUE    # Draw lines connecting each sample to its group centroid
BETA_SHOW_ELLIPSES <- FALSE   # Draw a confidence ellipse around each group
BETA_ELLIPSE_LEVEL <- 0.95    # Confidence level for ellipses

# --- 6b. Optional second grouping variable ---
# Set BETA_VARIABLE2_ENABLED <- TRUE and change the column name to use it
# in beta-diversity plots (as point shape).
BETA_VARIABLE2_ENABLED  <- TRUE
BETA_VARIABLE2_COLUMN   <- "MONTH"
BETA_VARIABLE2_AS_SHAPE <- TRUE   # Use Variable_2 as point shape in PCoA plots

# --- 6b-2. Optional permutation strata (blocking variable) ---
# Restricts PERMANOVA/betadisper permutations to occur only 
# WITHIN levels of this column instead of freely across the whole dataset
# This only configures WHICH column is available to use as strata, but is
# NOT applied automatically. Set `strata = TRUE` on individual models in
# BETA_PERMANOVA_MODELS (SECTION 6e) that should actually use it.
BETA_STRATA_ENABLED <- FALSE
BETA_STRATA_COLUMN  <- "COLONY"

# --- 6c. Distance methods ---
# source = "rel_abundance"  uses rarefied relative abundances
# source = "clean_counts"   uses cleaned, non-rarefied counts
# if FILE_TREE is NULL, wunifrac/unifrac are automatically dropped at run time
BETA_DISTANCE_METHODS <- list(
  list(name = "bray",     label = "Bray-Curtis",       source = "rel_abundance"),
  list(name = "aitchison",label = "Aitchison",         source = "clean_counts", pseudocount = 1),
  list(name = "wunifrac", label = "Weighted UniFrac",  source = "rel_abundance"),
  list(name = "unifrac",  label = "Unweighted UniFrac",source = "rel_abundance")
)

# --- 6d. PERMANOVA variables ---
VAR_1_PERMANOVA <- "GROUP"
VAR_2_PERMANOVA <- "MONTH"
VAR_3_PERMANOVA <- "FRAGMENT_STATUS"
VAR_4_PERMANOVA <- "SPECIES"

# --- 6e. PERMANOVA models ---
# `pairwise_var` (optional) selects which single variable is used for the
# pairwise post-hoc comparisons; defaults to VAR_1_PERMANOVA if omitted.
# `strata` Applies to both the PERMANOVA and its betadisper dispersion test.
BETA_PERMANOVA_MODELS <- list(
  species_only = list(
    enabled          = FALSE,
    rhs              = VAR_1_PERMANOVA,
    by               = NULL,
    dispersion_group = VAR_1_PERMANOVA,
    pairwise         = TRUE,
    pairwise_var     = VAR_1_PERMANOVA
  ),
  group_plus_species = list(
    enabled          = FALSE,
    rhs              = paste(VAR_2_PERMANOVA, "+", VAR_1_PERMANOVA),
    by               = "terms",
    dispersion_group = paste0("interaction(", VAR_2_PERMANOVA, ", ", VAR_1_PERMANOVA, ")"),
    pairwise         = FALSE
  ),
  month_x_species = list(
    enabled          = FALSE,
    rhs              = paste(VAR_2_PERMANOVA, "*", VAR_1_PERMANOVA),
    by               = "terms",
    dispersion_group = paste0("interaction(", VAR_2_PERMANOVA, ", ", VAR_1_PERMANOVA, ")"),
    pairwise         = FALSE
  ),
  group_plus_species_margin = list(
    enabled          = TRUE,
    rhs              = paste(VAR_1_PERMANOVA, "+", VAR_4_PERMANOVA),
    by               = "margin",
    dispersion_group = paste0("interaction(", VAR_1_PERMANOVA, ", ", VAR_4_PERMANOVA, ")"),
    pairwise         = FALSE
  ),
  lesion_vs_nolesion_paired = list(
    enabled          = FALSE,
    rhs              = VAR_3_PERMANOVA,
    by               = NULL,
    dispersion_group = VAR_3_PERMANOVA,
    pairwise         = FALSE,
    strata           = TRUE
  )
)

# --- 6f. Sub-feature toggles ---
BETA_PLOT_PCOA   <- TRUE   # Generate PCoA plots
BETA_PERMANOVA   <- TRUE   # Run PERMANOVA + betadisper
BETA_SAVE_MATRIX <- FALSE   # Export distance matrices as CSV
BETA_SAVE_COORDS <- FALSE   # Export PCoA coordinates as CSV

# ==============================================================================
# SECTION 7 — TAXONOMIC COMPOSITION
# ==============================================================================

TAXONOMY_ENABLED      <- TRUE   # Master switch for the entire taxonomy section

TAXONOMY_USE_RAREFIED <- TRUE
# TRUE  → taxonomy / relative abundance is computed on rarefied counts
# FALSE → taxonomy / relative abundance is computed on cleaned, non-rarefied counts

# --- 7a. Settings ---
TAXONOMY_RANK          <- "Family"
TAXONOMY_LABEL_COLUMN  <- "Family"
TAXONOMY_TOP_N         <- 19
TAXONOMY_SAMPLE_VAR    <- SAMPLE_ID_COLUMN   # x-axis label
TAXONOMY_PLOT_GROUP_VARS <- c(SAMPLE_ID_COLUMN) # TAXONOMY_SAMPLE_VAR must be included here

# Ordered list of paletteer palettes cycled across subsets.
TAXONOMY_PALETTES <- c(
  "ggthemes::Tableau_20",        # subset 1
  "ggsci::default_igv",          # subset 2
  "pals::alphabet",              # subset 3
  "Polychrome::palette36"        # subset 4
)

# Each entry uses its registered order from VARIABLE_LEVELS (SECTION 2) if
# one exists, otherwise it sorts alphabetically. Set to c() / NULL to just
# sort TAXONOMY_SAMPLE_VAR alphabetically.
TAXONOMY_SAMPLE_ORDER_VARS <- c("GROUP", "MONTH")

# --- 7b. Sub-feature toggles ---
TAXONOMY_PLOT     <- TRUE   # Generate stacked bar plots
TAXONOMY_SAVE_CSV <- TRUE   # Export all-taxa and top-taxa tables

# --- 7c — Canonical Relative Abundance Export 
# Independent of the stacked-bar taxonomy workflow above (SECTION 7). 
# Computed via the relative abundance helper file (see "Libraries" SECTION)
RELABUND_EXPORTS_ENABLED <- FALSE
RELABUND_EXPORT_SUBDIR   <- "relative_abundance"   # under OUTPUT_DIR/_global/

# ==============================================================================
# SECTION 8 — COLOR PALETTES & POINT SHAPES
# ==============================================================================

# Colors are applied consistently across ALL plots (alpha, beta, taxonomy).
# Set any palette to NULL to fall back to ggplot2 defaults for that variable.

PALETTE_GROUP <- c(
  "Acropora_Healthy_NoLesion"     = "#02A088",
  "Acropora_Diseased_NoLesion"    = "#F39B80",
  "Acropora_Diseased_Lesion"      = "#00AAD4",
  "Pocillopora_Healthy_NoLesion"  = "#8D5FD3",
  "Pocillopora_Diseased_NoLesion" = "#C87137",
  "Pocillopora_Diseased_Lesion"   = "#FFD42A"
)

PALETTE_MONTH <- c(
  "March" = "#1B7A8C",   
  "July"  = "#E07B39"    
)

# One fixed color per bacterial Family, for the Family-level taxonomic
# composition. Shared across every taxonomy subset  plot so the same family always gets the same color
PALETTE_FAMILY <- c(
  "Enterococcaceae"     = "#AC560E",
  "Flavobacteriaceae"   = "#E15759",
  "Cyclobacteriaceae"   = "#8FBC8F",
  "Oscillospiraceae"    = "#AFECDA",
  "Listeriaceae"        = "#4E79A7",
  "Simkaniaceae"        = "#A0CBE8",
  "Amoebophilaceae"     = "#66CDAA",
  "Sphingomonadaceae"   = "#B0C4DE",
  "Rikenellaceae"       = "#F28E2B",
  "Rhizobiaceae"        = "#8CD17D",
  "Myxococcaceae"       = "#B6992D",
  "Bacteroidaceae"      = "#F1CE63",
  "Staphylococcaceae"   = "#FFBE7D",
  "Pseudomonadaceae"    = "#59A14F",
  "Blattabacteriaceae"  = "#499894",
  "Muribaculaceae"      = "#86BCB6",
  "Arcobacteraceae"     = "#FF9D9A",
  "Lachnospiraceae"     = "#79706E",
  "Bacillaceae"         = "#FABFD2",
  "Lactobacillaceae"    = "#BAB0AC",
  "Paracoccaceae"       = "#D4A6C8",
  "Enterobacteriaceae"  = "#D37295",
  "Unclassified"        = "#B07AA1",
  "Endozoicomonadaceae" = "#9D7660",
  "Others"              = "#D7B5A6"
)

# Which palette above applies to each metadata column
PALETTE_BY_VARIABLE <- list(
  "GROUP"  = "group",
  "MONTH"  = "month",
  "Family" = "family"
)

# --- 8b. Point shapes ---
# A given value always renders with the same point shape in all plots.
# Set to NULL to fall back to ggplot2 default shapes.
SHAPE_MONTH <- c(
  "March" = 16,   # filled circle
  "July"  = 17    # filled triangle
)

# Which shape vector above applies to each metadata column.
SHAPE_BY_VARIABLE <- list(
  "MONTH" = "month"
)


# ==============================================================================
# RUN THE PIPELINE
# ==============================================================================

# Assemble flat config objects consumed by the pipeline functions
config <- list(
  project_name    = PROJECT_NAME,
  data_dir        = DATA_DIR,
  output_dir      = OUTPUT_DIR,
  files = list(
    taxonomy = FILE_TAXONOMY,
    metadata = FILE_METADATA,
    counts   = FILE_COUNTS,
    tree     = FILE_TREE
  ),
  sample_id_column = SAMPLE_ID_COLUMN,
  group_var        = GROUP_VAR,
  group_from       = GROUP_FROM,
  group_sep        = GROUP_SEP,
  group_levels     = GROUP_LEVELS,
  variable_levels  = VARIABLE_LEVELS,
  sample_group_label_column = SAMPLE_GROUP_LABEL_COLUMN,
  sample_group_label_sep    = SAMPLE_GROUP_LABEL_SEP,
  rarefaction = list(
    enabled     = RAREFACTION_ENABLED,
    sample_size = RAREFACTION_SAMPLE_SIZE,
    seed        = RAREFACTION_SEED,
    diagnostics = RAREFACTION_DIAGNOSTICS,
    curve_step  = RAREFACTION_CURVE_STEP
  ),
  subsets = SUBSETS,
  alpha = list(
    enabled         = ALPHA_ENABLED,
    use_rarefied    = ALPHA_USE_RAREFIED,
    metrics         = ALPHA_METRICS,
    plot_metrics    = ALPHA_PLOT_METRICS,
    y_limits        = ALPHA_Y_LIMITS,
    test_group_var  = ALPHA_GROUP_VAR,
    sig_threshold   = ALPHA_SIG_THRESHOLD,
    pairwise_adjust = ALPHA_PAIRWISE_ADJUST,
    plot            = ALPHA_PLOT,
    stats           = ALPHA_STATS,
    annotate        = ALPHA_ANNOTATE,
    show_points     = ALPHA_SHOW_POINTS,
    all_pvalues     = ALPHA_ALL_PVALUES,
    Variable_2 = list(
      enabled   = ALPHA_VARIABLE2_ENABLED,
      column    = ALPHA_VARIABLE2_COLUMN,
      as_facet  = ALPHA_VARIABLE2_AS_FACET,
      as_shape  = ALPHA_VARIABLE2_AS_SHAPE
    )
  ),
  beta = list(
    enabled            = BETA_ENABLED,
    use_rarefied       = BETA_USE_RAREFIED,
    ordination_method  = BETA_ORDINATION_METHOD,
    group_var          = BETA_GROUP_VAR,
    sample_label_var   = BETA_SAMPLE_LABEL_VAR,
    permutations       = BETA_PERMUTATIONS,
    seed               = BETA_SEED,
    pairwise_p_adjust  = BETA_PAIRWISE_P_ADJUST,
    show_segments      = BETA_SHOW_SEGMENTS,
    show_ellipses      = BETA_SHOW_ELLIPSES,
    ellipse_level      = BETA_ELLIPSE_LEVEL,
    Variable_2 = list(
      enabled    = BETA_VARIABLE2_ENABLED,
      column     = BETA_VARIABLE2_COLUMN,
      use_shape  = BETA_VARIABLE2_AS_SHAPE
    ),
    strata = list(
      enabled = BETA_STRATA_ENABLED,
      column  = BETA_STRATA_COLUMN
    ),
    distance_methods   = BETA_DISTANCE_METHODS,
    permanova_models   = BETA_PERMANOVA_MODELS,
    plot_pcoa          = BETA_PLOT_PCOA,
    run_permanova      = BETA_PERMANOVA,
    save_matrix        = BETA_SAVE_MATRIX,
    save_coords        = BETA_SAVE_COORDS
  ),
  taxonomy = list(
    enabled          = TAXONOMY_ENABLED,
    use_rarefied     = TAXONOMY_USE_RAREFIED,
    taxrank          = TAXONOMY_RANK,
    label_column     = TAXONOMY_LABEL_COLUMN,
    top_n            = TAXONOMY_TOP_N,
    sample_var       = TAXONOMY_SAMPLE_VAR,
    plot_group_vars  = TAXONOMY_PLOT_GROUP_VARS,
    sample_order_vars = TAXONOMY_SAMPLE_ORDER_VARS,
    palettes         = TAXONOMY_PALETTES,
    plot             = TAXONOMY_PLOT,
    save_csv         = TAXONOMY_SAVE_CSV
  ),
  relabund = list(
    enabled = RELABUND_EXPORTS_ENABLED,
    subdir  = RELABUND_EXPORT_SUBDIR
  ),
  palettes = list(
    group  = PALETTE_GROUP,
    month  = PALETTE_MONTH,
    family = PALETTE_FAMILY
  ),
  palette_by_variable = PALETTE_BY_VARIABLE,
  shapes = list(
    month = SHAPE_MONTH
  ),
  shape_by_variable = SHAPE_BY_VARIABLE
)


#### Utility functions ####

path_in <- function(..., cfg = config) file.path(cfg$data_dir, ...)
path_out <- function(..., cfg = config) file.path(cfg$output_dir, ...)
result_path <- function(subset_name, analysis_type, ..., cfg = config) {
  file.path(cfg$output_dir, subset_name, analysis_type, ...)
}

dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  invisible(path)
}

as_tibble_table <- function(x) {
  as.data.frame(x) %>%
    tibble::rownames_to_column("rowname")
}

capture_to_txt <- function(expr, file) {
  dir_create(dirname(file))
  capture.output(expr, file = file)
  invisible(file)
}

phyloseq_to_sample_matrix <- function(ps) {
  mat <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) mat <- t(mat)
  mat
}

drop_empty_taxa <- function(ps) {
  prune_taxa(taxa_sums(ps) > 0, ps)
}

filter_samples_by_rules <- function(ps, include = list(), exclude = list()) {
  meta <- data.frame(sample_data(ps), check.names = FALSE)
  keep <- rep(TRUE, nrow(meta))

  for (column in names(include)) {
    if (!column %in% colnames(meta)) stop("Missing metadata column: ", column, call. = FALSE)
    values <- include[[column]]
    keep <- keep & meta[[column]] %in% values
  }

  for (column in names(exclude)) {
    if (!column %in% colnames(meta)) stop("Missing metadata column: ", column, call. = FALSE)
    values <- exclude[[column]]
    keep <- keep & !(meta[[column]] %in% values)
  }

  ps %>%
    prune_samples(rownames(meta)[keep], .) %>%
    drop_empty_taxa()
}

save_csv <- function(x, file) {
  dir_create(dirname(file))
  x <- as.data.frame(x)
  if (!is.null(rownames(x)) && !all(rownames(x) == as.character(seq_len(nrow(x))))) {
    x <- tibble::rownames_to_column(x, "rowname")
  }
  readr::write_csv(x, file)
  invisible(file)
}

# Exports figures as PDF (vector format). Width and height are in inches.
save_plot <- function(plot, file, width = 8, height = 6) {
  # Ensure the extension is .pdf regardless of what the caller passes
  file <- sub("\\.(png|jpg|jpeg|tiff|svg)$", ".pdf", file, ignore.case = TRUE)
  if (!grepl("\\.pdf$", file, ignore.case = TRUE)) file <- paste0(file, ".pdf")
  dir_create(dirname(file))
  ggplot2::ggsave(file, plot = plot, width = width, height = height, device = "pdf")
  invisible(file)
}

safe_shapiro <- function(x) {
  x <- stats::na.omit(x)
  if (length(x) < 3 || length(x) > 5000) return(NULL)
  shapiro.test(x)
}

# Returns a ggplot2 discrete color/fill scale from a named palette vector,
# or NULL (ggplot2 default) when the palette is NULL or doesn't cover the data.
make_color_scale <- function(palette_vec, aesthetic = c("fill", "color"), drop = FALSE) {
  aesthetic <- match.arg(aesthetic)
  if (is.null(palette_vec)) return(NULL)
  if (aesthetic == "fill")
    ggplot2::scale_fill_manual(values = palette_vec, drop = drop, na.value = "grey70")
  else
    ggplot2::scale_color_manual(values = palette_vec, drop = drop, na.value = "grey70")
}

# Looks up the palette registered for a given metadata column (SECTION 8 — PALETTE_BY_VARIABLE). 
# Returns NULL (ggplot2 default colors) if the column has no registered palette
palette_for_variable <- function(cfg, var_name) {
  key <- cfg$palette_by_variable[[var_name]] %||% NULL
  if (is.null(key)) return(NULL)
  cfg$palettes[[key]] %||% NULL
}

# Same lookup as palette_for_variable(), but for the point-shape vectors
shape_for_variable <- function(cfg, var_name) {
  key <- cfg$shape_by_variable[[var_name]] %||% NULL
  if (is.null(key)) return(NULL)
  cfg$shapes[[key]] %||% NULL
}

# Returns a ggplot2 discrete shape scale from a named shape vector, or NULL
# (ggplot2 default shapes) when the vector is NULL or doesn't cover the data.
make_shape_scale <- function(shape_vec, drop = FALSE) {
  if (is.null(shape_vec)) return(NULL)
  ggplot2::scale_shape_manual(values = shape_vec, drop = drop, na.value = 1)
}

# TRUE when a tree file has been configured (SECTION 1 — FILE_TREE).
has_file_tree <- function(cfg) {
  tree_file <- cfg$files$tree
  !is.null(tree_file) && nzchar(tree_file)
}

# UniFrac-family distance methods require a phylogenetic tree.
distance_requires_tree <- function(method_name) {
  method_name %in% c("unifrac", "wunifrac")
}

#### Import, metadata, and phyloseq setup ####

load_microbiome_inputs <- function(cfg) {
  taxonomy <- read.csv(path_in(cfg$files$taxonomy, cfg = cfg), header = TRUE, row.names = 1, sep = "\t")
  metadata <- read.csv(path_in(cfg$files$metadata, cfg = cfg), header = TRUE, row.names = 1)
  counts   <- read.csv(path_in(cfg$files$counts,   cfg = cfg), header = TRUE, row.names = 1, sep = "\t")

  tree <- NULL
  if (has_file_tree(cfg)) {
    tree <- phyloseq::read_tree(path_in(cfg$files$tree, cfg = cfg))
    tree <- phytools::midpoint.root(tree)
  }

  list(
    taxonomy = as.matrix(taxonomy),
    metadata = data.frame(metadata, check.names = FALSE),
    counts   = as.matrix(counts),
    tree     = tree
  )
}

prepare_metadata <- function(metadata, cfg) {

  metadata[[cfg$sample_id_column]] <- rownames(metadata)

  if (length(cfg$group_from) > 0) {
    missing_group_cols <- setdiff(cfg$group_from, colnames(metadata))
    if (length(missing_group_cols) > 0) {
      stop("Missing metadata columns for group_from: ",
           paste(missing_group_cols, collapse = ", "), call. = FALSE)
    }
    metadata[[cfg$group_var]] <- do.call(
      paste, c(metadata[cfg$group_from], sep = cfg$group_sep)
    )
  }

  # Turn every column registered in VARIABLE_LEVELS (SECTION 2) into an ordered factor
  for (var_name in names(cfg$variable_levels)) {
    levels_vec <- cfg$variable_levels[[var_name]]
    if (is.null(levels_vec) || !var_name %in% colnames(metadata)) next
    metadata[[var_name]] <- ordered(metadata[[var_name]], levels = levels_vec)
  }

  if (!is.null(cfg$sample_group_label_column) && nzchar(cfg$sample_group_label_column)) {
    metadata[[cfg$sample_group_label_column]] <- paste(
      metadata[[cfg$sample_id_column]],
      metadata[[cfg$group_var]],
      sep = cfg$sample_group_label_sep %||% "_"
    )
  }

  metadata
}

build_phyloseq_object <- function(inputs, cfg) {
  metadata <- prepare_metadata(inputs$metadata, cfg)
  components <- list(
    otu_table(inputs$counts, taxa_are_rows = TRUE),
    tax_table(inputs$taxonomy),
    sample_data(metadata)
  )
  if (!is.null(inputs$tree)) components <- c(components, list(phy_tree(inputs$tree)))
  do.call(phyloseq, components)
}

clean_phyloseq_object <- function(ps) {
  ps %>%
    subset_taxa(Kingdom == "Bacteria" & Family != "Mitochondria" & Order != "Chloroplast") %>%
    prune_taxa(taxa_sums(.) > 1, .)
}

run_rarefaction_diagnostics <- function(ps, cfg) {
  diag_dir <- result_path("_global", "rarefaction", cfg = cfg)
  dir_create(diag_dir)

  # Reads per sample
  reads_per_sample <- data.frame(sum = sample_sums(ps))
  save_csv(reads_per_sample, file.path(diag_dir, "reads_per_sample.csv"))

  # Rarefaction curve: ASVs discovered vs. sequencing depth, one line per
  # sample, with a vertical line at the chosen rarefaction threshold.
  otu_mat    <- phyloseq_to_sample_matrix(ps)
  curve_file <- file.path(diag_dir, "rarefaction_curve.pdf")
  grDevices::pdf(curve_file, width = 10, height = 8)
  vegan::rarecurve(otu_mat, step = cfg$rarefaction$curve_step %||% 10000,
                    label = TRUE, col = "blue",
                    xlab = "Sequencing depth", ylab = "ASVs")
  graphics::abline(v = cfg$rarefaction$sample_size, col = "red", lty = 2)
  grDevices::dev.off()

  invisible(list(reads_per_sample = reads_per_sample, curve_file = curve_file))
}

rarefy_phyloseq_object <- function(ps, cfg) {
  if (!cfg$rarefaction$enabled) return(ps)
  set.seed(cfg$rarefaction$seed)
  rarefy_even_depth(
    ps,
    sample.size = cfg$rarefaction$sample_size,
    replace     = FALSE,
    rngseed     = cfg$rarefaction$seed
  )
}

make_relative_abundance <- function(ps, multiplier = 100) {
  transform_sample_counts(ps, function(x) multiplier * x / sum(x))
}

make_subset_collection <- function(ps, subset_specs) {
  purrr::imap(subset_specs, function(spec, name) {
    filter_samples_by_rules(
      ps,
      include = spec$include %||% list(),
      exclude = spec$exclude %||% list()
    )
  })
}

#### Alpha diversity ####

calculate_alpha_diversity <- function(ps, metrics, group_var) {
  richness <- estimate_richness(ps) %>%
    tibble::rownames_to_column(".sample_row") %>%
    left_join(
      data.frame(sample_data(ps), check.names = FALSE) %>%
        tibble::rownames_to_column(".sample_row"),
      by = ".sample_row"
    )

  if ("Evenness" %in% metrics) {
    richness <- richness %>%
      mutate(Evenness = ifelse(Observed > 1, Shannon / log(Observed), NA_real_))
  }

  summary <- richness %>%
    group_by(.data[[group_var]]) %>%
    summarise(
      across(
        any_of(metrics),
        list(
          mean   = ~ mean(.x,   na.rm = TRUE),
          median = ~ median(.x, na.rm = TRUE),
          sd     = ~ sd(.x,     na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    )

  list(richness = richness, summary = summary)
}

plot_alpha_metric <- function(richness, metric, group_var, title = metric,
                              y_limits = NULL, palette = NULL,
                              sig_pairs = NULL,
                              show_points = FALSE,
                              all_pvalues = FALSE,
                              facet_var = NULL,
                              shape_var = NULL,
                              shape_palette = NULL) {
  p <- richness %>%
    ggplot(aes(x = .data[[group_var]], y = .data[[metric]], fill = .data[[group_var]])) +
    geom_boxplot(outlier.shape = NA) +
    theme_bw() +
    theme(
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      axis.title        = element_text(size = 16),
      axis.text.y       = element_text(size = 14),
      axis.text.x       = element_text(angle = 45, hjust = 1),
      legend.position   = "right",
      text              = element_text(size = 14)
    ) +
    labs(x = NULL, y = metric, fill = group_var, title = title)

  scale <- make_color_scale(palette, "fill")
  if (!is.null(scale)) p <- p + scale

  # Facet by a second grouping variable (e.g. MONTH) — one panel per level
  if (!is.null(facet_var)) {
    p <- p + facet_wrap(facets = vars(.data[[facet_var]]))
  }

  # Overlay individual sample dots (jittered to avoid overlap)
  if (isTRUE(show_points)) {
    if (!is.null(shape_var)) {
      p <- p + geom_jitter(
        aes(color = .data[[group_var]], shape = .data[[shape_var]]),
        width = 0.15, size = 1.8, alpha = 0.7,
        show.legend = c(color = FALSE, shape = TRUE)
      )
      shape_scale <- make_shape_scale(shape_palette)
      if (!is.null(shape_scale)) p <- p + shape_scale
    } else {
      p <- p + geom_jitter(
        aes(color = .data[[group_var]]),
        width = 0.15, size = 1.8, alpha = 0.7, show.legend = FALSE
      )
    }
    color_scale <- make_color_scale(palette, "color")
    if (!is.null(color_scale)) p <- p + color_scale
  }

  # Add brackets: either all pairwise p-values or significant-only (* / ** / ***)
  if (!is.null(sig_pairs) && nrow(sig_pairs) > 0) {
    if (isTRUE(all_pvalues)) {
      bracket_data <- sig_pairs %>%
        mutate(label = formatC(p_adjusted, format = "f", digits = 3)) %>%
        select(group1, group2, label)
    } else {
      bracket_data <- sig_pairs %>%
        filter(!is.na(stars)) %>%
        select(group1, group2, label = stars)
    }

    if (nrow(bracket_data) > 0) {
      # Extend y-axis so brackets don't overlap the data
      y_max   <- max(richness[[metric]], na.rm = TRUE)
      y_range <- diff(range(richness[[metric]], na.rm = TRUE))
      y_step  <- y_range * 0.08
      y_start <- y_max + y_step

      bracket_data <- bracket_data %>%
        mutate(y.position = y_start + (seq_len(nrow(.)) - 1) * y_step)

      p <- p +
        ggpubr::stat_pvalue_manual(
          bracket_data,
          label          = "label",
          y.position     = "y.position",
          tip.length     = 0.01,
          bracket.size   = 0.4,
          size           = 4
        )
    }
  }

  if (!is.null(y_limits)) p <- p + coord_cartesian(ylim = y_limits)
  p
}

run_alpha_tests <- function(richness, metrics, group_var,
                            sig_threshold  = 0.05,
                            pairwise_adjust = "BH") {

  # ── helper: p-value → asterisk label (or NA to suppress) ──────────────────
  p_to_stars <- function(p) {
    dplyr::case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < sig_threshold ~ "*",
      TRUE      ~ NA_character_
    )
  }

  purrr::map(metrics, function(metric) {

    values <- richness[[metric]]
    groups <- richness[[group_var]]

    # ── 1. Shapiro-Wilk per group ──────────────────────────────────────────
    sw_results <- richness %>%
      group_by(.data[[group_var]]) %>%
      summarise(
        n       = sum(!is.na(.data[[metric]])),
        SW_W    = {
          x <- stats::na.omit(.data[[metric]])
          if (length(x) >= 3 && length(x) <= 5000) unname(shapiro.test(x)$statistic) else NA_real_
        },
        SW_p    = {
          x <- stats::na.omit(.data[[metric]])
          if (length(x) >= 3 && length(x) <= 5000) shapiro.test(x)$p.value else NA_real_
        },
        .groups = "drop"
      ) %>%
      mutate(normal = SW_p >= sig_threshold | is.na(SW_p))

    normality_passed <- all(sw_results$normal, na.rm = TRUE)

    # ── 2. Levene's test (robust to non-normality) ─────────────────────────
    levene_result <- tryCatch(
      car::leveneTest(values ~ factor(groups), center = median),
      error = function(e) NULL
    )
    levene_p      <- if (!is.null(levene_result)) levene_result$`Pr(>F)`[1] else NA_real_
    homogeneity_passed <- !is.na(levene_p) && levene_p >= sig_threshold

    # ── 3. Decision: parametric or non-parametric ─────────────────────────
    use_parametric <- normality_passed && homogeneity_passed

    method_note <- dplyr::case_when(
      use_parametric                         ~ "ANOVA + Tukey (normality OK, homogeneity OK)",
      normality_passed & !homogeneity_passed ~ "Kruskal-Wallis + Dunn (normality OK but heteroscedastic)",
      !normality_passed                      ~ "Kruskal-Wallis + Dunn (normality failed)"
    )

    # ── 4a. ANOVA + Tukey (parametric path) ───────────────────────────────
    anova_result <- tukey_result <- tukey_pairs <- NULL
    if (use_parametric) {
      formula     <- stats::as.formula(paste(metric, "~", group_var))
      anova_model <- aov(formula, data = richness)
      anova_result <- summary(anova_model)
      anova_p      <- anova_result[[1]][["Pr(>F)"]][1]

      if (!is.na(anova_p) && anova_p < sig_threshold) {
        tukey_result <- tryCatch(TukeyHSD(anova_model, group_var), error = function(e) NULL)
        if (!is.null(tukey_result)) {
          tukey_pairs <- as.data.frame(tukey_result[[group_var]]) %>%
            tibble::rownames_to_column("comparison") %>%
            rename(p_adjusted = `p adj`) %>%
            mutate(
              group1 = sub("^(.+)-(.+)$", "\\1", comparison),
              group2 = sub("^(.+)-(.+)$", "\\2", comparison),
              stars  = p_to_stars(p_adjusted)
            )
        }
      }
    }

    # ── 4b. Kruskal-Wallis + Dunn (non-parametric path) ───────────────────
    kruskal_result <- dunn_pairs <- NULL
    if (!use_parametric) {
      formula        <- stats::as.formula(paste(metric, "~", group_var))
      kruskal_result <- kruskal.test(formula, data = richness)

      if (!is.na(kruskal_result$p.value) && kruskal_result$p.value < sig_threshold) {
        dunn_raw <- dunn.test::dunn.test(
          values, groups,
          method = pairwise_adjust,
          list   = FALSE, table = FALSE, kw = FALSE, label = FALSE
        )
        dunn_pairs <- tibble(
          comparison = dunn_raw$comparisons,
          Z          = dunn_raw$Z,
          p_adjusted = dunn_raw$P.adjusted
        ) %>%
          mutate(
            group1 = trimws(sub("^(.+) - (.+)$", "\\1", comparison)),
            group2 = trimws(sub("^(.+) - (.+)$", "\\2", comparison)),
            stars  = p_to_stars(p_adjusted)
          )
      }
    }

    # ── 5. Collect significant pairs for plot annotation ─────────────────
    sig_pairs <- if (use_parametric) tukey_pairs else dunn_pairs
    sig_pairs_plot <- if (!is.null(sig_pairs)) {
      sig_pairs %>% filter(!is.na(stars))
    } else {
      NULL
    }

    list(
      metric              = metric,
      use_parametric      = use_parametric,
      method_note         = method_note,
      shapiro_per_group   = sw_results,
      levene              = levene_result,
      levene_p            = levene_p,
      anova               = anova_result,
      tukey               = tukey_result,
      tukey_pairs         = tukey_pairs,
      kruskal             = kruskal_result,
      dunn_pairs          = dunn_pairs,
      sig_pairs           = sig_pairs,
      sig_pairs_plot      = sig_pairs_plot
    )
  }) %>%
    rlang::set_names(metrics)
}

save_alpha_tests <- function(tests, alpha_dir) {
  # ── per-metric summary table ───────────────────────────────────────────────
  summary_rows <- purrr::map_dfr(tests, function(t) {
    omnibus_p <- if (t$use_parametric) {
      tryCatch(t$anova[[1]][["Pr(>F)"]][1], error = function(e) NA_real_)
    } else {
      t$kruskal$p.value %||% NA_real_
    }
    tibble(
      metric          = t$metric,
      method          = t$method_note,
      levene_p        = t$levene_p,
      omnibus_p       = omnibus_p,
      n_sig_pairs     = if (!is.null(t$sig_pairs_plot)) nrow(t$sig_pairs_plot) else 0L
    )
  })
  save_csv(summary_rows, file.path(alpha_dir, "alpha_test_summary.csv"))

  # ── per-metric detail files ────────────────────────────────────────────────
  purrr::walk(tests, function(t) {
    m <- t$metric

    # Shapiro per group
    save_csv(t$shapiro_per_group,
             file.path(alpha_dir, paste0(m, "_shapiro_per_group.csv")))

    # Levene
    if (!is.null(t$levene)) {
      capture_to_txt(print(t$levene),
                     file.path(alpha_dir, paste0(m, "_levene.txt")))
    }

    # Omnibus (whichever was run)
    if (t$use_parametric && !is.null(t$anova)) {
      capture_to_txt(print(t$anova),
                     file.path(alpha_dir, paste0(m, "_anova.txt")))
    }
    if (!t$use_parametric && !is.null(t$kruskal)) {
      capture_to_txt(print(t$kruskal),
                     file.path(alpha_dir, paste0(m, "_kruskal.txt")))
    }

    # Post-hoc pairwise (whichever was run, if omnibus was significant)
    if (t$use_parametric && !is.null(t$tukey_pairs)) {
      save_csv(t$tukey_pairs,
               file.path(alpha_dir, paste0(m, "_tukey_pairs.csv")))
    }
    if (!t$use_parametric && !is.null(t$dunn_pairs)) {
      save_csv(t$dunn_pairs,
               file.path(alpha_dir, paste0(m, "_dunn_pairs.csv")))
    }
  })
}

# TRUE when a Variable_2 column has been configured for alpha diversity
# (SECTION 5d) — mirrors beta_variable_2_column() for the beta section.
alpha_variable_2_column <- function(cfg) {
  variable_2 <- cfg$alpha$Variable_2 %||% list(enabled = FALSE)
  if (!isTRUE(variable_2$enabled)) return(NULL)
  variable_2$column %||% NULL
}

run_alpha_workflow <- function(subsets, cfg) {
  variable_2_column <- alpha_variable_2_column(cfg)
  facet_var <- if (!is.null(variable_2_column) && isTRUE(cfg$alpha$Variable_2$as_facet)) variable_2_column else NULL
  shape_var <- if (!is.null(variable_2_column) && isTRUE(cfg$alpha$Variable_2$as_shape)) variable_2_column else NULL

  purrr::imap(subsets, function(ps, subset_name) {
    message("Alpha diversity: ", subset_name)
    alpha_dir <- result_path(subset_name, "alpha_diversity", cfg = cfg)
    dir_create(alpha_dir)

    alpha <- calculate_alpha_diversity(ps, cfg$alpha$metrics, cfg$alpha$test_group_var)

    save_csv(alpha$richness, file.path(alpha_dir, "alpha_richness.csv"))
    save_csv(alpha$summary,  file.path(alpha_dir, "alpha_summary.csv"))

    # Run assumption-aware statistical tests first (needed for annotations)
    tests <- NULL
    if (isTRUE(cfg$alpha$stats)) {
      tests <- run_alpha_tests(
        alpha$richness,
        metrics         = cfg$alpha$metrics,
        group_var       = cfg$alpha$test_group_var,
        sig_threshold   = cfg$alpha$sig_threshold   %||% 0.05,
        pairwise_adjust = cfg$alpha$pairwise_adjust %||% "BH"
      )
      save_alpha_tests(tests, alpha_dir)
    }

    if (isTRUE(cfg$alpha$plot)) {
      purrr::walk(cfg$alpha$plot_metrics, function(metric) {
        title <- if (metric == "Observed") "ASV Count" else paste(metric, "Index")

        # Pull pairs for annotation (if enabled and stats were run)
        sig_pairs <- NULL
        if (isTRUE(cfg$alpha$annotate) && !is.null(tests) && metric %in% names(tests)) {
          if (isTRUE(cfg$alpha$all_pvalues)) {
            # Use ALL pairwise pairs so every comparison gets a p-value label
            sig_pairs <- tests[[metric]]$sig_pairs
          } else {
            sig_pairs <- tests[[metric]]$sig_pairs_plot
          }
        }

        p <- plot_alpha_metric(
          alpha$richness,
          metric      = metric,
          group_var   = cfg$alpha$test_group_var,
          title       = paste(subset_name, title),
          y_limits    = cfg$alpha$y_limits[[metric]],
          palette     = palette_for_variable(cfg, cfg$alpha$test_group_var),
          sig_pairs   = sig_pairs,
          show_points = isTRUE(cfg$alpha$show_points),
          all_pvalues = isTRUE(cfg$alpha$all_pvalues),
          facet_var     = facet_var,
          shape_var     = shape_var,
          shape_palette = if (!is.null(shape_var)) shape_for_variable(cfg, shape_var) else NULL
        )
        save_plot(p, file.path(alpha_dir, paste0(metric, ".pdf")), width = 8, height = 6)
      })
    }

    c(alpha, list(tests = tests))
  })
}

#### Beta diversity ####

make_distance_matrix <- function(ps, method, pseudocount = 1) {
  if (method == "aitchison") {
    mat             <- phyloseq_to_sample_matrix(ps)
    distance_matrix <- vegan::vegdist(mat, method = "aitchison", pseudocount = pseudocount)
  } else {
    distance_matrix <- phyloseq::distance(ps, method = method)
  }
  stopifnot(identical(sample_names(ps), rownames(as.matrix(distance_matrix))))
  distance_matrix
}

beta_variable_2_column <- function(cfg) {
  variable_2 <- cfg$beta$Variable_2 %||% list(enabled = FALSE)
  if (!isTRUE(variable_2$enabled)) return(NULL)
  variable_2$column %||% NULL
}

# The strata COLUMN NAME for a PERMANOVA model (SECTION 6b-2 and SECTION 6e), 
# or NULL if strata aren't configured/enabled or this particular model doesn't request them.
beta_strata_column <- function(cfg, model_cfg) {
  if (!isTRUE(model_cfg$strata)) return(NULL)
  strata_cfg <- cfg$beta$strata %||% list(enabled = FALSE)
  if (!isTRUE(strata_cfg$enabled)) return(NULL)
  strata_cfg$column %||% NULL
}

# Builds a permute::how() control object restricting permutations to occur only WITHIN levels of `strata_vec`
make_permutation_control <- function(n_permutations, strata_vec = NULL) {
  if (is.null(strata_vec)) return(n_permutations)
  permute::how(nperm = n_permutations, blocks = strata_vec)
}

validate_beta_variables <- function(metadata, cfg) {
  required           <- cfg$beta$group_var
  variable_2_column  <- beta_variable_2_column(cfg)
  if (!is.null(variable_2_column)) required <- c(required, variable_2_column)

  missing <- setdiff(required, colnames(metadata))
  if (length(missing) > 0)
    stop("Missing beta-diversity metadata column(s): ", paste(missing, collapse = ", "), call. = FALSE)

  invisible(TRUE)
}

make_beta_metadata <- function(ps, cfg) {
  metadata <- data.frame(sample_data(ps), check.names = FALSE)
  validate_beta_variables(metadata, cfg)
  metadata$Group <- metadata[[cfg$beta$group_var]]

  variable_2_column <- beta_variable_2_column(cfg)
  if (!is.null(variable_2_column)) metadata$Variable_2 <- metadata[[variable_2_column]]

  metadata
}

run_pcoa <- function(ps, distance_matrix, method = "PCoA", group_var = "Group", sample_label_var = "Sample_ID") {
  ordination <- ordinate(ps, method = method, distance = distance_matrix)
  coords     <- as.data.frame(ordination$vectors)
  colnames(coords)[1:2] <- c("PC1", "PC2")

  meta <- data.frame(sample_data(ps), check.names = FALSE) %>%
    tibble::rownames_to_column(".sample_row")

  coords <- coords %>%
    tibble::rownames_to_column(".sample_row") %>%
    left_join(meta, by = ".sample_row")

  if (any(is.na(coords[[group_var]])))
    stop("Metadata did not join correctly to ordination coordinates.", call. = FALSE)

  centroids <- coords %>%
    group_by(.data[[group_var]]) %>%
    summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop")

  coords_with_centroids <- coords %>%
    left_join(centroids, by = group_var, suffix = c("", "_centroid"))

  eig_vals    <- ordination$values$Eigenvalues
  var_explained <- round(100 * eig_vals / sum(eig_vals), 1)

  list(
    ordination            = ordination,
    coordinates           = coords,
    centroids             = centroids,
    coords_with_centroids = coords_with_centroids,
    var_explained         = var_explained,
    group_var             = group_var,
    sample_label_var      = sample_label_var
  )
}

plot_pcoa <- function(pcoa_result, title, distance_label, shape_var = NULL,
                      palette = NULL, shape_palette = NULL,
                      show_segments = TRUE, show_ellipses = FALSE,
                      ellipse_level = 0.95) {
  group_var <- pcoa_result$group_var

  p <- ggplot(pcoa_result$coordinates, aes(x = PC1, y = PC2, color = .data[[group_var]]))

  if (!is.null(shape_var) && shape_var %in% colnames(pcoa_result$coordinates))
    p <- p + aes(shape = .data[[shape_var]])

  if (isTRUE(show_segments)) {
    p <- p +
      geom_segment(
        data = pcoa_result$coords_with_centroids,
        aes(x = PC1_centroid, y = PC2_centroid, xend = PC1, yend = PC2,
            color = .data[[group_var]]),
        linewidth = 0.8, alpha = 0.5
      )
  }

  if (isTRUE(show_ellipses)) {
    p <- p +
      stat_ellipse(
        aes(group = .data[[group_var]], color = .data[[group_var]]),
        level = ellipse_level, linewidth = 0.6
      )
  }

  p <- p +
    geom_point(size = 3.5) +
    geom_point(
      data       = pcoa_result$centroids,
      aes(x = PC1, y = PC2, color = .data[[group_var]]),
      shape      = 4, size = 4, stroke = 2,
      inherit.aes = FALSE
    ) +
    labs(
      x     = paste0("PC1 (", pcoa_result$var_explained[1], "%)"),
      y     = paste0("PC2 (", pcoa_result$var_explained[2], "%)"),
      color = group_var,
      title = paste(title, "-", distance_label)
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position   = "right",
      panel.grid.major  = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor  = element_line(color = "grey95", linewidth = 0.2),
      plot.title        = element_text(face = "bold")
    )

  scale <- make_color_scale(palette, "color")
  if (!is.null(scale)) p <- p + scale

  shape_scale <- make_shape_scale(shape_palette)
  if (!is.null(shape_scale)) p <- p + shape_scale
  p
}

permanova_formula <- function(rhs) stats::as.formula(paste("distance_matrix ~", rhs))

# Extracts the metadata column names referenced in a formula-like string
# (rhs, dispersion_group, or a single variable name), so they can be
# validated against the metadata up front.
extract_formula_vars <- function(expr_str) {
  if (is.null(expr_str)) return(character(0))
  all.vars(parse(text = expr_str)[[1]])
}

# Evaluates any rhs/dispersion_group expression (a bare column name, or an
# expression like "interaction(Treatment_mother, Group)") directly against
# the metadata data frame.
dispersion_groups <- function(metadata, dispersion_group) {
  dispersion_group <- dispersion_group %||% "Group"
  tryCatch(
    eval(parse(text = dispersion_group)[[1]], envir = metadata),
    error = function(e) {
      stop("Could not evaluate dispersion_group '", dispersion_group, "': ",
           conditionMessage(e), call. = FALSE)
    }
  )
}

pairwise_permanova <- function(distance_matrix, metadata, group_var,
                               permutations = 999, seed = 69, p_adjust = "BH",
                               strata_column = NULL) {
  groups <- unique(metadata[[group_var]])
  groups <- groups[!is.na(groups)]
  pairs  <- combn(as.character(groups), 2, simplify = FALSE)

  set.seed(seed)
  result <- purrr::map_dfr(pairs, function(pair) {
    sample_ids <- rownames(metadata)[metadata[[group_var]] %in% pair]
    sub_dist   <- as.dist(as.matrix(distance_matrix)[sample_ids, sample_ids])
    sub_meta   <- metadata[sample_ids, , drop = FALSE]
    sub_meta$.group <- droplevels(factor(sub_meta[[group_var]]))
    # Re-derive strata from sub_meta (not the full metadata) — the restricted
    # permutation scheme must be built from whichever samples actually
    # survive this pair's subset.
    sub_strata   <- if (!is.null(strata_column)) sub_meta[[strata_column]] else NULL
    perm_control <- make_permutation_control(permutations, sub_strata)
    fit <- vegan::adonis2(sub_dist ~ .group, data = sub_meta, permutations = perm_control)
    tibble(
      pair      = paste(pair, collapse = " vs "),
      Df        = fit$Df[1],
      SumOfSqs  = fit$SumOfSqs[1],
      R2        = fit$R2[1],
      F_model   = fit$F[1],
      p_value   = fit$`Pr(>F)`[1]
    )
  })

  result %>% mutate(p_adjusted = p.adjust(p_value, method = p_adjust))
}

run_permanova_model <- function(distance_matrix, metadata, model_cfg,
                                permutations = 999, seed = 69, p_adjust = "BH",
                                strata_column = NULL) {
  stopifnot(identical(rownames(metadata), rownames(as.matrix(distance_matrix))))

  strata_vec   <- if (!is.null(strata_column)) metadata[[strata_column]] else NULL
  perm_control <- make_permutation_control(permutations, strata_vec)

  set.seed(seed)
  dispersion      <- betadisper(distance_matrix, dispersion_groups(metadata, model_cfg$dispersion_group %||% "Group"))
  dispersion_test <- permutest(dispersion, permutations = perm_control)

  set.seed(seed)
  permanova <- vegan::adonis2(
    permanova_formula(model_cfg$rhs %||% "Group"),
    data         = metadata,
    permutations = perm_control,
    by           = model_cfg$by %||% NULL
  )

  pairwise <- NULL
  if (isTRUE(model_cfg$pairwise)) {
    pairwise <- pairwise_permanova(
      distance_matrix = distance_matrix,
      metadata        = metadata,
      group_var       = model_cfg$pairwise_var %||% "Group",
      permutations    = permutations,
      seed            = seed,
      p_adjust        = p_adjust,
      strata_column   = strata_column
    )
  }

  list(
    dispersion      = dispersion,
    dispersion_test = dispersion_test,
    permanova       = permanova,
    pairwise        = pairwise
  )
}

run_permanova <- function(ps, distance_matrix, cfg) {
  metadata       <- make_beta_metadata(ps, cfg)
  models         <- cfg$beta$permanova_models %||% list(group_only = list(enabled = TRUE, rhs = "Group"))
  enabled_models <- purrr::keep(models, ~ isTRUE(.x$enabled))

  if (length(enabled_models) == 0) return(list())

  # Fail fast with a clear message if any enabled model references a
  # metadata column that doesn't exist (e.g. a typo in a VAR_*_PERMANOVA
  # definition, or a column that was renamed).
  required_vars <- unique(unlist(purrr::map(enabled_models, function(m) {
    c(extract_formula_vars(m$rhs %||% "Group"),
      extract_formula_vars(m$dispersion_group %||% "Group"),
      extract_formula_vars(m$pairwise_var %||% NULL),
      beta_strata_column(cfg, m))
  })))
  missing_vars <- setdiff(required_vars, colnames(metadata))
  if (length(missing_vars) > 0) {
    stop("PERMANOVA model(s) reference metadata column(s) not found: ",
         paste(missing_vars, collapse = ", "),
         ". Check the VAR_*_PERMANOVA definitions, BETA_PERMANOVA_MODELS, ",
         "and BETA_STRATA_COLUMN.",
         call. = FALSE)
  }

  purrr::imap(enabled_models, function(model_cfg, model_name) {
    run_permanova_model(
      distance_matrix = distance_matrix,
      metadata        = metadata,
      model_cfg       = model_cfg,
      permutations    = cfg$beta$permutations,
      seed            = cfg$beta$seed,
      p_adjust        = cfg$beta$pairwise_p_adjust,
      strata_column   = beta_strata_column(cfg, model_cfg)
    )
  })
}

save_permanova_results <- function(permanova_results, prefix) {
  purrr::iwalk(permanova_results, function(model_result, model_name) {
    model_prefix <- paste(prefix, model_name, sep = "_")
    capture_to_txt(model_result$dispersion_test, paste0(model_prefix, "_betadisper.txt"))
    capture_to_txt(model_result$permanova,        paste0(model_prefix, "_permanova.txt"))
    if (!is.null(model_result$pairwise))
      save_csv(model_result$pairwise, paste0(model_prefix, "_pairwise_permanova.csv"))

    if (model_name == "group_only") {
      capture_to_txt(model_result$dispersion_test, paste0(prefix, "_betadisper.txt"))
      capture_to_txt(model_result$permanova,        paste0(prefix, "_permanova.txt"))
      if (!is.null(model_result$pairwise))
        save_csv(model_result$pairwise, paste0(prefix, "_pairwise_permanova.csv"))
    }
  })
}

run_beta_workflow <- function(source_subsets, cfg) {
  purrr::imap(source_subsets, function(subset_sources, subset_name) {
    message("Beta diversity: ", subset_name)
    beta_dir <- result_path(subset_name, "beta_diversity", cfg = cfg)
    dir_create(beta_dir)

    purrr::map(cfg$beta$distance_methods, function(method_cfg) {
      method_name  <- method_cfg$name
      method_label <- method_cfg$label %||% method_name
      source_name  <- method_cfg$source %||% "rel_abundance"
      ps           <- subset_sources[[source_name]]

      if (is.null(ps))
        stop("Missing beta source '", source_name, "' for subset ", subset_name, call. = FALSE)

      distance_matrix <- make_distance_matrix(ps, method = method_name,
                                              pseudocount = method_cfg$pseudocount %||% 1)

      prefix <- file.path(beta_dir, method_name)

      if (isTRUE(cfg$beta$save_matrix))
        save_csv(as.matrix(distance_matrix), paste0(prefix, "_distance_matrix.csv"))

      pcoa <- NULL
      p    <- NULL
      if (isTRUE(cfg$beta$plot_pcoa)) {
        pcoa <- run_pcoa(
          ps, distance_matrix  = distance_matrix,
          method          = cfg$beta$ordination_method,
          group_var       = cfg$beta$group_var,
          sample_label_var = cfg$beta$sample_label_var
        )

        shape_var         <- NULL
        variable_2_column <- beta_variable_2_column(cfg)
        if (isTRUE((cfg$beta$Variable_2 %||% list())$use_shape) && !is.null(variable_2_column))
          shape_var <- variable_2_column

        p <- plot_pcoa(pcoa, cfg$subsets[[subset_name]]$title %||% subset_name,
                       method_label, shape_var = shape_var,
                       palette       = palette_for_variable(cfg, cfg$beta$group_var),
                       shape_palette = if (!is.null(shape_var)) shape_for_variable(cfg, shape_var) else NULL,
                       show_segments = isTRUE(cfg$beta$show_segments %||% TRUE),
                       show_ellipses = isTRUE(cfg$beta$show_ellipses %||% FALSE),
                       ellipse_level = cfg$beta$ellipse_level %||% 0.95)
        save_plot(p, paste0(prefix, "_pcoa.pdf"), width = 8, height = 6)

        if (isTRUE(cfg$beta$save_coords)) {
          save_csv(pcoa$coordinates, paste0(prefix, "_pcoa_coordinates.csv"))
          save_csv(pcoa$centroids,   paste0(prefix, "_pcoa_centroids.csv"))
        }
      }

      permanova_results <- NULL
      if (isTRUE(cfg$beta$run_permanova)) {
        permanova_results <- run_permanova(ps, distance_matrix, cfg)
        save_permanova_results(permanova_results, prefix)
      }

      list(
        method    = method_name,
        distance  = distance_matrix,
        pcoa      = pcoa,
        plot      = p,
        permanova = permanova_results
      )
    }) %>%
      rlang::set_names(purrr::map_chr(cfg$beta$distance_methods, "name"))
  })
}


#### Taxonomic composition ####

make_taxonomy_table <- function(ps, cfg) {
  taxrank <- cfg$taxonomy$taxrank

  taxonomy <- ps %>%
    tax_glom(taxrank = taxrank) %>%
    transform_sample_counts(function(x) x / sum(x)) %>%
    psmelt() %>%
    filter(Abundance != 0) %>%
    mutate(across(where(is.factor), as.character))

  if (taxrank == "Genus" && all(c("Genus", "Family") %in% colnames(taxonomy))) {
    taxonomy <- taxonomy %>%
      mutate(
        Genus_Family = ifelse(
          is.na(Genus) | Genus %in% c("Unclassified", "endosymbionts", ""),
          paste(Genus, Family, sep = "_"),
          Genus
        )
      )
  }

  taxonomy
}

taxonomy_sample_order <- function(ps, sample_var, order_vars, variable_levels = list()) {
  meta       <- data.frame(sample_data(ps), check.names = FALSE)
  order_vars <- intersect(order_vars %||% character(0), colnames(meta))

  rank_cols <- purrr::map(order_vars, function(v) {
    levels_vec <- variable_levels[[v]]
    if (!is.null(levels_vec)) match(meta[[v]], levels_vec) else as.integer(factor(meta[[v]]))
  })

  ord <- do.call(order, c(rank_cols, list(meta[[sample_var]])))
  unique(meta[[sample_var]][ord])
}

top_taxa_by_sample <- function(taxonomy, label_column, sample_var, plot_group_vars, top_n = 19) {
  top_names <- taxonomy %>%
    group_by(.data[[label_column]]) %>%
    summarise(Total_Abundance = sum(Abundance), .groups = "drop") %>%
    arrange(desc(Total_Abundance)) %>%
    slice_head(n = top_n) %>%
    mutate(is_top = TRUE)

  taxonomy_grouped <- taxonomy %>%
    left_join(top_names %>% select(all_of(label_column), is_top), by = label_column) %>%
    mutate(.taxon_label = ifelse(is.na(is_top), "Others", .data[[label_column]]))

  taxonomy_grouped[[label_column]] <- taxonomy_grouped$.taxon_label

  taxonomy_grouped %>%
    select(-.taxon_label, -is_top) %>%
    group_by(.data[[sample_var]]) %>%
    mutate(totalSum = sum(Abundance)) %>%
    ungroup() %>%
    group_by(across(any_of(c(plot_group_vars, label_column)))) %>%
    summarise(
      New_TaxonAbundance = sum(Abundance),
      totalSum           = unique(totalSum),
      Relative_Abundance = New_TaxonAbundance / unique(totalSum),
      .groups            = "drop"
    ) %>%
    unique()
}

plot_taxonomic_composition <- function(taxonomy_top, label_column, sample_var,
                                       palette, group_palette = NULL, title,
                                       sample_levels = NULL, manual_palette = NULL) {
  abundance_rank <- taxonomy_top %>%
    group_by(.data[[label_column]]) %>%
    summarise(total_abundance = sum(Relative_Abundance), .groups = "drop") %>%
    arrange(desc(total_abundance))

  bottom_labels <- intersect(c("Others", "Unclassified_Unclassified"), abundance_rank[[label_column]])
  rest_labels   <- setdiff(abundance_rank[[label_column]], bottom_labels)   # already desc-sorted
  bottom_to_top <- c(bottom_labels, rest_labels)

  taxonomy_top[[label_column]] <- factor(taxonomy_top[[label_column]], levels = rev(bottom_to_top))

  # Order the x-axis (see taxonomy_sample_order() — TAXONOMY_SAMPLE_ORDER_VARS in SECTION 7). 
  if (!is.null(sample_levels)) {
    taxonomy_top[[sample_var]] <- factor(taxonomy_top[[sample_var]], levels = sample_levels)
  }

  p <- ggplot(taxonomy_top) +
    geom_col(aes(x = .data[[sample_var]], y = Relative_Abundance, fill = .data[[label_column]]),
             color = "black") +
    ylab("Proportion of Community") +
    xlab(NULL) +
    theme_linedraw() +
    theme(
      axis.text.y     = element_text(size = 20, color = "black"),
      axis.title.y    = element_text(size = 18, color = "black"),
      axis.text.x     = element_text(size = 15, angle = 60, vjust = 1, hjust = 1, color = "black"),
      legend.text     = element_text(size = 10),
      legend.position = "right",
      strip.text      = element_text(size = 10, face = "bold"),
      legend.title    = element_text(face = "bold", size = 14),
      plot.title      = element_text(face = "bold")
    ) +
    labs(fill = label_column, title = title) +
    guides(fill = guide_legend(ncol = 1, byrow = TRUE))

  fill_scale <- make_color_scale(manual_palette, "fill", drop = TRUE)
  if (!is.null(fill_scale)) {
    p + fill_scale
  } else {
    p + paletteer::scale_fill_paletteer_d(palette)
  }
}

run_taxonomy_workflow <- function(subsets, cfg) {
  subset_names <- names(subsets)
  palettes     <- cfg$taxonomy$palettes %||% "ggthemes::Tableau_20"

  purrr::imap(subsets, function(ps, subset_name) {
    message("Taxonomic composition: ", subset_name)
    taxonomy_dir <- result_path(subset_name, "taxonomic_composition", cfg = cfg)
    dir_create(taxonomy_dir)

    # Pick palette by subset position, cycling if fewer palettes than subsets
    subset_index   <- which(subset_names == subset_name)
    palette_index  <- ((subset_index - 1L) %% length(palettes)) + 1L
    subset_palette <- palettes[[palette_index]]

    taxonomy     <- make_taxonomy_table(ps, cfg)
    label_column <- cfg$taxonomy$label_column
    if (!label_column %in% colnames(taxonomy)) label_column <- cfg$taxonomy$taxrank

    taxonomy_top <- top_taxa_by_sample(
      taxonomy,
      label_column    = label_column,
      sample_var      = cfg$taxonomy$sample_var,
      plot_group_vars = cfg$taxonomy$plot_group_vars,
      top_n           = cfg$taxonomy$top_n
    )

    if (isTRUE(cfg$taxonomy$save_csv)) {
      save_csv(taxonomy,     file.path(taxonomy_dir, paste0(cfg$taxonomy$taxrank, "_all.csv")))
      save_csv(taxonomy_top, file.path(taxonomy_dir, paste0(cfg$taxonomy$taxrank, "_top_taxa.csv")))
    }

    p <- NULL
    if (isTRUE(cfg$taxonomy$plot)) {
      sample_levels <- taxonomy_sample_order(
        ps,
        sample_var      = cfg$taxonomy$sample_var,
        order_vars      = cfg$taxonomy$sample_order_vars,
        variable_levels = cfg$variable_levels
      )

      p <- plot_taxonomic_composition(
        taxonomy_top,
        label_column   = label_column,
        sample_var     = cfg$taxonomy$sample_var,
        palette        = subset_palette,
        group_palette  = cfg$palettes$group,
        title          = paste(cfg$subsets[[subset_name]]$title %||% subset_name, "-", cfg$taxonomy$taxrank),
        sample_levels  = sample_levels,
        manual_palette = palette_for_variable(cfg, label_column)
      )
      save_plot(p, file.path(taxonomy_dir, paste0(cfg$taxonomy$taxrank, "_top_taxa.pdf")),
                width = 12, height = 7)
    }

    list(all = taxonomy, top = taxonomy_top, plot = p)
  })
}


#### Canonical relative abundance exports (SECTION 7c) ####

# Per-sample and group mean/min/max relative abundance at Family and Genus
# level, produced through the shared helpers (shushah_relabund_helpers.R). 
#`ps_source` is the rarefied or the cleaned non-rarefied phyloseq object, chosen from cfg$taxonomy$use_rarefied.
run_relabundance_exports <- function(ps_source, cfg, use_rarefied) {
  out_dir <- result_path("_global", cfg$relabund$subdir %||% "relative_abundance", cfg = cfg)
  dir_create(out_dir)

  meta <- data.frame(sample_data(ps_source), check.names = FALSE)
  # coral genus x health-status (AH/DU/DL), months pooled
  meta$STATUS_ABBR <- factor(
    dplyr::recode(
      paste(meta$COLONY_STATUS, meta$FRAGMENT_STATUS, sep = "_"),
      "Healthy_NoLesion"  = "AH",
      "Diseased_NoLesion" = "DU",
      "Diseased_Lesion"   = "DL"
    ),
    levels = c("AH", "DU", "DL")
  )
  sample_meta <- data.frame(
    SampleID    = rownames(meta),
    GENUS       = as.character(meta$GENUS),
    STATUS_ABBR = meta$STATUS_ABBR,
    check.names = FALSE
  )
  group_cols <- c("GENUS", "STATUS_ABBR")

  for (level in c("Family", "Genus")) {
    long  <- shushah_relabund_long(ps_source, level = level)
    summ  <- shushah_group_summary(long, sample_meta, group_cols = group_cols)
    top19 <- shushah_topN_union_summary(long, sample_meta, group_cols = group_cols,
                                        subset_col = "GENUS", top_n = 19)
    readr::write_csv(long,  file.path(out_dir, paste0(level, "_relabundance_persample.csv")))
    readr::write_csv(summ,  file.path(out_dir, paste0(level, "_relabundance_group_summary.csv")))
    readr::write_csv(top19, file.path(out_dir, paste0(level, "_relabundance_top19_group_summary.csv")))
    message("  ", level, ": ", nrow(long), " per-sample rows, ", nrow(summ),
            " group rows, ", length(attr(top19, "named_taxa")), " named + Others (top-19 union)")
  }

  shushah_write_relab_methods(
    file.path(out_dir, "RELAB_METHODS.txt"),
    use_rarefied = use_rarefied,
    group_desc   = paste(
      "coral genus x health-status (AH = Healthy/NoLesion, DU = Diseased/NoLesion,",
      "DL = Diseased/Lesion); March + July pooled"
    ),
    extra_lines  = c(
      "Source script: ShushahTotal_Microbiome_Composition_PipelineTest3.R (SECTION 7c)",
      "This is the project's source for relative-abundance values."
    )
  )
  message("Canonical relative abundance exports written to: ", out_dir)
  invisible(out_dir)
}


#### Main run ####

run_microbiome_composition_pipeline <- function(cfg = config) {
  dir_create(cfg$output_dir)

  inputs       <- load_microbiome_inputs(cfg)

  # No tree loaded → drop any UniFrac distance methods rather than erroring
  # deep inside phyloseq::distance().
  if (is.null(inputs$tree)) {
    tree_dependent <- purrr::keep(cfg$beta$distance_methods, ~ distance_requires_tree(.x$name))
    if (length(tree_dependent) > 0) {
      message("No phylogenetic tree loaded — skipping UniFrac distance method(s): ",
              paste(purrr::map_chr(tree_dependent, "name"), collapse = ", "))
      cfg$beta$distance_methods <- purrr::keep(cfg$beta$distance_methods, ~ !distance_requires_tree(.x$name))
    }
  }

  raw_ps       <- build_phyloseq_object(inputs, cfg)
  clean_ps     <- clean_phyloseq_object(raw_ps)

  if (isTRUE(cfg$rarefaction$diagnostics)) {
    message("Rarefaction diagnostics: reads-per-sample + rarefaction curve")
    run_rarefaction_diagnostics(clean_ps, cfg)
  }

  rarefied_ps  <- rarefy_phyloseq_object(clean_ps, cfg)

  # SECTION 7c — canonical relative abundance exports
  if (isTRUE(cfg$relabund$enabled)) {
    use_rar <- isTRUE(cfg$taxonomy$use_rarefied)
    message("Canonical relative abundance exports: ",
            if (use_rar) "rarefied" else "clean non-rarefied", " counts")
    run_relabundance_exports(if (use_rar) rarefied_ps else clean_ps, cfg,
                             use_rarefied = use_rar)
  }

  # Relative abundance versions — one from rarefied, one from clean counts
  rel_ab_rarefied_ps <- make_relative_abundance(rarefied_ps,  multiplier = 100)
  rel_ab_clean_ps    <- make_relative_abundance(clean_ps,     multiplier = 100)

  # Subset collections for every combination needed
  clean_count_subsets    <- make_subset_collection(clean_ps,           cfg$subsets)
  rarefied_count_subsets <- make_subset_collection(rarefied_ps,        cfg$subsets)
  rel_ab_rarefied_subsets <- make_subset_collection(rel_ab_rarefied_ps, cfg$subsets)
  rel_ab_clean_subsets    <- make_subset_collection(rel_ab_clean_ps,    cfg$subsets)

  beta_rel_ab_subsets <- if (isTRUE(cfg$beta$use_rarefied)) {
    rel_ab_rarefied_subsets
  } else {
    rel_ab_clean_subsets
  }

  beta_sources <- purrr::imap(cfg$subsets, function(spec, subset_name) {
    list(
      clean_counts    = clean_count_subsets[[subset_name]],
      rarefied_counts = rarefied_count_subsets[[subset_name]],
      rel_abundance   = beta_rel_ab_subsets[[subset_name]]
    )
  })

  sample_summary <- tibble(
    object   = c("raw", "clean", "rarefied"),
    ntaxa    = c(ntaxa(raw_ps),    ntaxa(clean_ps),    ntaxa(rarefied_ps)),
    nsamples = c(nsamples(raw_ps), nsamples(clean_ps), nsamples(rarefied_ps))
  )
  save_csv(sample_summary,
           result_path("_global", "pipeline", "phyloseq_object_summary.csv", cfg = cfg))

  results <- list(
    phyloseq = list(raw = raw_ps, clean = clean_ps, rarefied = rarefied_ps,
                    rel_abundance_rarefied = rel_ab_rarefied_ps,
                    rel_abundance_clean    = rel_ab_clean_ps),
    subsets  = list(clean_counts          = clean_count_subsets,
                    rarefied_counts       = rarefied_count_subsets,
                    rel_abundance_rarefied = rel_ab_rarefied_subsets,
                    rel_abundance_clean    = rel_ab_clean_subsets),
    alpha    = NULL,
    beta     = NULL,
    taxonomy = NULL
  )

  if (isTRUE(cfg$alpha$enabled)) {
    alpha_subsets <- if (isTRUE(cfg$alpha$use_rarefied)) {
      rarefied_count_subsets
    } else {
      clean_count_subsets
    }
    message("Alpha diversity input: ",
            if (isTRUE(cfg$alpha$use_rarefied)) "rarefied counts" else "clean non-rarefied counts")
    results$alpha <- run_alpha_workflow(alpha_subsets, cfg)
  }

  if (isTRUE(cfg$beta$enabled)) {
    message("Beta diversity rel_abundance source: ",
            if (isTRUE(cfg$beta$use_rarefied)) "rarefied" else "clean non-rarefied")
    results$beta <- run_beta_workflow(beta_sources, cfg)
  }

  if (isTRUE(cfg$taxonomy$enabled)) {
    taxonomy_subsets <- if (isTRUE(cfg$taxonomy$use_rarefied)) {
      rarefied_count_subsets
    } else {
      clean_count_subsets
    }
    message("Taxonomy input: ",
            if (isTRUE(cfg$taxonomy$use_rarefied)) "rarefied counts" else "clean non-rarefied counts")
    results$taxonomy <- run_taxonomy_workflow(taxonomy_subsets, cfg)
  }

  dir_create(result_path("_global", "pipeline", cfg = cfg))
  saveRDS(results, result_path("_global", "pipeline",
                               "microbiome_composition_results.rds", cfg = cfg))
  results
}

results <- run_microbiome_composition_pipeline(config)
