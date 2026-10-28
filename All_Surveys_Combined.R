#### Libraries ####

library(ggplot2)
library(dplyr)
library(scales)
library(tidyr)
library(readr)
library(forcats)
library(ggsci)
library(paletteer)

#### Setting everything up ####

# Set theme for graphs
theme_set(theme_bw())

# Set directory
setwd(
  '/Users/casartm/Desktop/KAUST PhD/Projects/Coral Disease - Shushah/Disease Surveys/All Analyses Combined')

# Read files
CoralCover <- read.csv("CoralCover_Combined.csv")
colnames(CoralCover)

PIT <- read.csv("PIT_Combined.csv")
colnames(PIT)
unique(PIT$Classes)

Diseases <- read.csv("CoralDiseases_Combined.csv")
colnames(Diseases)

#### Add metadata and prepare file - Coral Cover ####

# Create a new column combining Site and Transect with an underscore
CoralCover$SiteTransect <-
  paste(CoralCover$Site, CoralCover$Transect, sep = "_")
head(CoralCover)

# Create a new column with the reef category information
CoralCover$Reef_category <-
  fct_collapse(
    CoralCover$Site,
    Restoration = c("Coral_Head" ,
                    "Pinnacle_B",
                    "Crest_E"),
    Control = c("Crest_A"),
    External_control = c("NW_Shushah")
  )

# Remove unusable sites from analysis
CoralCover <- CoralCover %>%
  filter(!(Site %in% c("Coral_Head")))

# Create a metadata file
CoralCover_metadata <- CoralCover %>%
  select(Site, Reef_category, Length) %>%
  distinct()

# Calculate total area covered per site
CoralCover_area <- CoralCover %>%
  select(Site, SiteTransect, Length) %>%
  distinct() %>%  # Removes duplicates caused by genus repetition
  group_by(Site) %>%
  summarise(SiteArea = sum(Length), .groups = "drop")

# Collapse genera
CoralCover$Genus <-
  fct_collapse(
    CoralCover$Genus,
    Millepora = c("Millepora", "Millepora*"),
    Porites = c("Porites", "Porites*")
  )

# Store the size classes (columns) names in a vector
size_columns <-
  c(
    "under5cm",
    "X5to10cm",
    "X10to20cm",
    "X20to40cm",
    "X40to80cm",
    "X80to160cm",
    "above160cm"
  )

#### Add metadata and prepare file - PIT ####

# Create a new column combining Site and Transect with an underscore
PIT$SiteTransect <- 
  paste(PIT$Site, PIT$Transect, sep = "_")
head(PIT)

# Remove unusable sites from alaysis
PIT <- PIT %>%
  filter(!(Site %in% c("Coral_Head", "Hadda_Beach")))

# Define classes colour
classes_palette <- c(
  CCA = "#f4a4c4",     # pastel pink
  DC  = "#aeeae0",     # aquamarine
  HC  = "#d2b48c",     # light brown (tan)
  MA  = "#b4e7b4",     # pastel green
  MI  = "#ffb347",     # pastel orange
  OC  = "#add8e6",     # light blue
  OT  = "#cfa5e2",     # pastel purple
  RO  = "#cccccc",     # pastel grey
  RU  = "#d5c6f0",     # pastel lilac
  SA  = "#fff799",     # pastel yellow
  SP  = "#ff9999",     # pastel red
  TU  = "#77dd77"      # pastel dark green
)

#### Add metadata and prepare file - Diseases ####

# Remove unusable sites from alaysis
Diseases <- Diseases %>%
  filter(!(Site %in% c("Haddah_Beach")))

# Create a metadata file
Diseases_metadata <- Diseases %>%
  select(Site, Transect_Area) %>%
  distinct()

# Convert to long format
Diseases_long <- Diseases %>%
  pivot_longer(cols = -c(Site, Transect_Area),
               names_to = "Genus_Disease",
               values_to = "Count")

# Create metadata table by separating Genus and Disease Type
Diseases_genera <- Diseases_long %>%
  select(Genus_Disease) %>%
  distinct() %>%  # Keep only unique entries
  separate(
    Genus_Disease,
    into = c("Genus", "Disease_Type"),
    sep = "_",
    remove = FALSE
  )

#### CoralCover - Calculate coral cover values ####

# Calculate the total number of corals per transect for each genus (row)
CoralCover$TotalGenus <- rowSums(CoralCover[, size_columns])
head(CoralCover)

# Calculate the coral density per square meter
CoralCover$CoralDensity <- CoralCover$TotalGenus / CoralCover$Length
head(CoralCover)
# According to ChatGPT, this measure is best suited to compare coral density across sites
# accounting for the differences in sampling effort (10 vs 20 m^2)

# Count all corals per transect and aggregate in a new dataset with a single row per transect
CoralCover_summary <- CoralCover %>%
  group_by(Site, SiteTransect) %>%
  summarize(
    TotalCorals = sum(TotalGenus),
    # Total coral count for the transect
    .groups = "drop"               # Prevent nested groupings
  )

# Add reef category
CoralCover_summary <- CoralCover_summary %>%
  left_join(CoralCover_metadata, by = c("Site"))

# Calculate the coral density at each transect
CoralCover_summary$CoralDensityTransect <-
  CoralCover_summary$TotalCorals / CoralCover_summary$Length
head(CoralCover_summary)
unique(CoralCover_summary)

# Calculate average density
CoralCover_AverageDensity_Site <- CoralCover_summary %>%
  group_by(Site) %>%
  summarise(MeanDensity = mean(CoralDensityTransect)) 

# Filter out genera with no observations
CoralCover_filtered <- CoralCover %>%
  filter(TotalGenus > 0)
unique(CoralCover_filtered$Genus)

# Obtain the overall proportion of each genera
CoralCover_GenusProportion <-
  CoralCover_filtered %>%
  group_by(Genus) %>%
  summarise(Total_byGenera = sum(TotalGenus)) %>%
  ungroup() %>%
  mutate(Genus_proportion = Total_byGenera / sum(Total_byGenera))

# Collapse genera with overall proportion < 1% under "Others"
CoralCover_GenusProportion <-
  CoralCover_GenusProportion %>%
  mutate(
    Genus_collapsed = ifelse(
      Genus_proportion < 0.01 |
        Genus == "Unknown",
      "Others",
      as.character(Genus)
    )
  )

# Left-join the names
CoralCover_filtered <- 
  CoralCover_filtered %>%
  left_join(CoralCover_GenusProportion %>%
              select(Genus, Genus_collapsed),
            by = "Genus")
unique(CoralCover_filtered$Genus_collapsed)

# Left-join the site area
CoralCover_filtered <-
  CoralCover_filtered %>%
  left_join(CoralCover_area,
            by = "Site"
  )
unique(CoralCover_filtered$SiteArea)

# AS IT IS JUST DESCRIPTIVE, I SUM UP ALL TRANSECTS AND DIVIDE BY THE TOTAL LENGHT
# INSTEAD OF DOING THE MEANS

# Summarize total counts by genus and site
CoralCover_RelAb <- CoralCover_filtered %>%
  group_by(Site, Genus_collapsed, SiteArea) %>%
  summarize(TotalCount = sum(TotalGenus), .groups = "drop") 

# Calculate the coral density for each genera
CoralCover_RelAb <- CoralCover_RelAb %>%
  group_by(Site, Genus_collapsed) %>%
  mutate(CoralDensity_byGenus = TotalCount / SiteArea)

# Normalize to get proportions
CoralCover_RelAb <- CoralCover_RelAb %>%
  group_by(Site) %>%
  mutate(Proportion = TotalCount / sum(TotalCount)) %>%
  ungroup()

# Calculate total density by sum of different genera
CoralCover_RelAb %>%
  group_by(Site) %>%
  summarise(TotalSiteDensity = sum(CoralDensity_byGenus))

#### PIT - Calculate substrate classes values ####

# Point-interval columns are named by their distance along the transect (e.g. "X0.25", "X0.5", ..., "X20")
point_columns <- grep("^X[0-9.]+$", names(PIT), value = TRUE)

# Calculate the total number of observations per transect for each class (row)
PIT$TotalClass <- rowSums(PIT[, point_columns])
head(PIT)

# Double-check the total of points per transect
PIT %>%
  group_by(SiteTransect) %>%
  summarise(TotalPoints = sum(TotalClass)) %>%
  print(n = Inf)

# Calculate the total transect length
PIT_TransectLength <- PIT %>%
  group_by(SiteTransect) %>%
  summarize(TransectLength = sum(TotalClass) * 0.25)

# Calculate the proportion and percentage of each class in a new file
PIT_summary <- PIT %>%
  group_by(Site) %>%
  mutate(TotalSiteCount = sum(TotalClass)) %>%
  ungroup() %>%
  group_by(Site, Classes) %>%
  summarize(
    Class_Proportion = sum(TotalClass) / unique(TotalSiteCount),
    Class_Percentage = Class_Proportion * 100,
    .groups = "drop"
  )

#### Diseases - Calculate FOC and Prevalence ####

# Calculate total disease count per site
Diseases_summary_site <- Diseases_long %>%
  group_by(Site) %>%
  summarize(TotalDisease_Site = sum(Count))

### Disease FOC ###

# Calculate disease occurrence
Diseases_FOC <- Diseases_long %>%
  filter(Count > 0) %>%               # Keep only rows where the disease is present
  distinct(Site, Genus_Disease) %>%   # One row per Site × Genus_Disease pair
  count(Genus_Disease, name = "Disease_occurrence")  # Count number of sites per Genus_Disease

# Left-join coral genus and disease type
Diseases_FOC <- Diseases_FOC %>%
  left_join(Diseases_genera, by = "Genus_Disease")

# Extract unique genera with diseases
Diseases_DiseasedGenera <- unique(Diseases_FOC$Genus)
Diseases_DiseasedGenera

# Filter genera from the coral cover dataset that were reported with diseases
CoralCover_DiseasedGenera_count <- CoralCover_filtered %>%
  filter(Genus %in% Diseases_DiseasedGenera) %>%
  distinct(Site, Genus)

# Count the number of sites with diseased genera present
CoralCover_DiseasedGenera_count <- CoralCover_DiseasedGenera_count %>%
  count(Genus, name = "Sites_with_genus")

# Merge to the disease occurrence dataset and calculate FOC
Diseases_FOC <- Diseases_FOC %>%
  left_join(CoralCover_DiseasedGenera_count, by = "Genus") %>%
  mutate(FOC = Disease_occurrence / Sites_with_genus)

### Disease Prevalence by Site ###

# Filter genus count data by retaining only genera with disease reported
CoralCover_TotalGenus_Diseased <- CoralCover_filtered %>%
  filter(Genus %in% Diseases_DiseasedGenera) %>%
  select(Genus, Site, SiteTransect, TotalGenus, SiteArea)

# Obtain the density of each genus at each site
CoralCover_SiteGenusDensity_Diseased <- CoralCover_TotalGenus_Diseased %>%
  group_by(Site, Genus, SiteArea) %>%
  summarize(Site_Count = sum(TotalGenus), .groups = "drop")

# Left-join to the long-format disease data and calculate prevalence by site
Diseases_prevalence_bySite <- Diseases_long %>%
  filter(Count > 0) %>%
  left_join(Diseases_genera, by = c("Genus_Disease")) %>%
  left_join(CoralCover_SiteGenusDensity_Diseased, by = c("Site", "Genus")) %>%
  mutate(Prevalence = ((Count / Transect_Area) / (Site_Count / SiteArea)) * 100)

## Calculate disease prevalence at each site, disregarding disease type

# Here, I need to divide the total disease count per site by the total number of corals observed

# Calculate the total coral count per site
CoralCover_SiteSummary <- CoralCover_summary %>%
  group_by(Site) %>%
  summarize(TotalCoralCount = sum(TotalCorals))

# Left-join to the disease summary by site and caculate prevalence
Diseases_prevalence_SitesOnly <- Diseases_summary_site %>%
  left_join(CoralCover_SiteSummary, by = c("Site")) %>%
  left_join(Diseases_metadata, by = c("Site")) %>%
  left_join(CoralCover_area, by = c("Site")) %>%
  mutate(Prevalence = ((TotalDisease_Site / Transect_Area) / (TotalCoralCount / SiteArea)) * 100)
  

### Diseases prevalence - overall ###

## 1. Here, I am first calculating the total transect area where the genera of interest were found
## 2. Then I am calculating the total overall count of these genera
## 3. Then, I am calculating the total area surveyed for diseases where the genera of interest were found
## 4. And finally, I am calculating the total overall count for each disease
## The overall disease prevalence will be Prevalence_overall = [(4/3) / (2/1)] * 100

# 1. Calculate the total transect area where the genus affected by diseases were found
CoralCover_TotalDiseaseArea <- CoralCover_SiteGenusDensity_Diseased %>%
  group_by(Genus) %>%
  summarize(Total_CoralCover_Area = sum(SiteArea))

# 2. Recover values of total count for the genera affected by diseases
Diseases_genera_TotalCount <- CoralCover_GenusProportion %>%
  filter(Genus %in% Diseases_DiseasedGenera) %>%
  select(Genus, Total_byGenera)

# 3. Calculate the total transect area where diseases were surveyed 
Diseases_TotalDiseaseArea <- CoralCover_SiteGenusDensity_Diseased %>%
  select(Site, Genus) %>%
  left_join(Diseases_metadata, by = "Site") %>%
  group_by(Genus) %>%
  summarize(Total_Disease_Area = sum(Transect_Area))

# Summarize by disease type and calculate overall prevalence
Diseases_prevalence_overall <- Diseases_prevalence_bySite %>%
  group_by(Genus_Disease) %>%
  summarise(Overall_Count_Disease = sum(Count),
            .groups = "drop",) %>%
  left_join(Diseases_genera, by ="Genus_Disease") %>%
  left_join(Diseases_genera_TotalCount, by = "Genus") %>%
  left_join(CoralCover_TotalDiseaseArea, by = "Genus") %>%
  left_join(Diseases_TotalDiseaseArea, by = "Genus") %>%
  mutate(
    Total_Disease_Area = Total_Disease_Area,
    Prevalence = ((Overall_Count_Disease / Total_Disease_Area) /
                    (Total_byGenera / Total_CoralCover_Area)
    ) * 100
  )

# Calculate disease count by genus (not by disease type)
Diseases_prevalence_overall_byGenus <- Diseases_prevalence_overall %>%
  group_by(Genus) %>%
  summarise(Overall_Count_Disease_byGenus = sum(Overall_Count_Disease),
            .groups = "drop")

# Calculate overall prevalence by genus (not by disease type)
Diseases_prevalence_overall_byGenus <- Diseases_prevalence_overall_byGenus %>%
  left_join(Diseases_genera_TotalCount, by = "Genus") %>%
  left_join(CoralCover_TotalDiseaseArea, by = "Genus") %>%
  left_join(Diseases_TotalDiseaseArea, by = "Genus") %>%
  mutate(Prevalence_byGenus = ((Overall_Count_Disease_byGenus / Total_Disease_Area) /
                                 (Total_byGenera / Total_CoralCover_Area)
  ) * 100)

#### CoralCover - Plot coral relative abundances ####

# Reorder ranking abundances from lower to higher
CoralCover_RelAb <- CoralCover_RelAb %>%
  mutate(Genus_collapsed = fct_reorder(Genus_collapsed, Proportion, sum))

# Get the desired factor levels from the other object
genus_levels <- levels(CoralCover_RelAb$Genus_collapsed)

# Plotting relative abundances by site
ggplot(CoralCover_RelAb, aes(x = Site, y = Proportion, fill = Genus_collapsed)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Coral Composition by Genus and Site",
       x = "Site",
       y = "Proportion of Coral Cover",
       fill = "Coral Genus") +
  theme_minimal() +
  scale_fill_d3(palette = "category20")

# Plotting total count by site
ggplot(CoralCover_RelAb, aes(x = Site, y = TotalCount, fill = Genus_collapsed)) +
  geom_bar(stat = "identity") +
  # scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Coral Composition by Genus and Site",
       x = "Site",
       y = "Coral Cover",
       fill = "Coral Genus") +
  theme_minimal() +
  scale_fill_d3(palette = "category20c")

# Plotting coral densities
ggplot(CoralCover_RelAb, aes(x = Site, y = CoralDensity_byGenus, fill = Genus_collapsed)) +
  geom_bar(stat = "identity") +
  # scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Coral Composition by Genus and Site",
       x = "Site",
       y = "Coral density",
       fill = "Coral Genus") +
  theme_minimal() +
  scale_fill_paletteer_d("ggthemes::Tableau_20")

#### CoralCover - Plot coral cover by size classes ####

# Calculate the total count per size class for each transect using pivot_longer
CoralCover_SizeClass <- CoralCover_filtered %>%
  select(Site, SiteTransect, Length, size_columns, Genus_collapsed) %>%
  pivot_longer(cols = size_columns, names_to = "SizeClass", values_to = "Count") %>%
  left_join(CoralCover_area, by = "Site") %>%
  group_by(Site, Genus_collapsed, SizeClass, SiteArea) %>%
  summarise(TotalCount = sum(Count), .groups = "drop") %>%
  mutate(Density = TotalCount / SiteArea)
head(CoralCover_SizeClass)

# Reorder genera by total relative abundance
CoralCover_SizeClass <- CoralCover_SizeClass %>%
  mutate(Genus_collapsed = factor(Genus_collapsed, levels = genus_levels))

# Reorder size classes 
CoralCover_SizeClass <- CoralCover_SizeClass %>%
  mutate(SizeClass = factor(SizeClass, levels = size_columns))

# Plot coral cover by size classes (facet by Site)
ggplot(CoralCover_SizeClass,
       aes(x = SizeClass, y = Density, fill = Genus_collapsed)) +
  geom_col(position = "stack") +
  facet_wrap( ~ Site) +
  labs(title = "Coral Density by Genus and Size Class",
       y = "Density (count per m²)",
       x = "Size Class") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_paletteer_d("ggthemes::Tableau_20")

# Calculate the size class proportion at each site
CoralCover_SizeClass_Proportions <- CoralCover_filtered %>%
  select(Site, SiteTransect, Length, size_columns, Genus_collapsed) %>%
  pivot_longer(cols = size_columns, names_to = "SizeClass", values_to = "Count") %>%
  group_by(Site, SizeClass) %>%
  summarise(SizeClass_Count = sum(Count), .groups = "drop") %>%
  group_by(Site) %>%
  mutate(SizeClass_Proportion = SizeClass_Count / sum(SizeClass_Count)) %>%
  ungroup()
head(CoralCover_SizeClass_Proportions)

# Reorder size classes
CoralCover_SizeClass_Proportions <-
  CoralCover_SizeClass_Proportions %>%
  mutate(SizeClass = factor(SizeClass, levels = rev(size_columns)))

# Plot coral cover by size classes (facet by Site)
ggplot(CoralCover_SizeClass_Proportions,
       aes(x = Site, y = SizeClass_Proportion, fill = SizeClass)) +
  geom_col(position = "stack") +
  labs(title = "Proportion of Size Classes",
       y = "Size class proportion",
       x = "Site") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set3")

#### Plot substrate composition####

# Re-order proportions
PIT_summary <- PIT_summary %>%
  mutate(Classes = fct_reorder(Classes, Class_Proportion, sum))
head(PIT_summary)

# Plot substrate composition by site proportion
ggplot(PIT_summary, aes(x = Site, y = Class_Proportion, fill = Classes)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Site Substrate Composition",
       x = "Site",
       y = "Substrate Proportion",
       fill = "Substrate Class") +
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = classes_palette)

#### Diseases - Plot metrics ####

# Create vector to order disease genus types
DiseaseGenus_type_levels <-
  c (
    "Alveopora_FocalBleaching",
    "Millepora_FocalBleaching",
    "Porites_FocalBleaching",
    
    "Acropora_GA",
    "Echinopora_GA",
    "Goniastrea_GA",
    "Montipora_GA",
    "Porites_GA",
    
    "Porites_MultifocalPinkSpot",
    
    "Acropora_WS",
    "Cyphastrea_WS",
    "Millepora_WS",
    "Montipora_WS",
    "Pocillopora_WS",
    "Porites_WS",
    "Stylophora_WS",
    
    "Goniastrea_WSandSEB",
    "Pocillopora_WSandSEB"
  )

#Re-order disease genus types
Diseases_prevalence_overall$Genus_Disease <-
  factor(Diseases_prevalence_overall$Genus_Disease,
         levels = DiseaseGenus_type_levels)
levels(Diseases_prevalence_overall$Genus_Disease)

# Plot overall disease prevalence by disease type
ggplot(
  Diseases_prevalence_overall,
  aes(x = Genus_Disease, y = Prevalence, fill = Disease_Type)
) +
  geom_col() +
  coord_cartesian(ylim = c(0, max(
    Diseases_prevalence_overall$Prevalence[Diseases_prevalence_overall$Genus_Disease != "Echinopora_GA"]
  ) * 1)) +
  labs(x = "Disease Type and Genus", y = "Prevalence") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")

# Plot overall disease prevalence by site
ggplot(
  Diseases_prevalence_SitesOnly,
  aes(x = Site, y = Prevalence, fill = Site)
) +
  geom_col() +
  # coord_cartesian(ylim = c(0, max(
  #   Diseases_prevalence_overall$Prevalence[Diseases_prevalence_overall$Genus_Disease != "Echinopora_GA"]
  # ) * 1)) +
  labs(x = "Site", y = "Disease Prevalence") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")

#### Print useful CSV files ####

# Print Diseases FOC
write_csv(
  Diseases_FOC %>% 
    select(Genus, 
           Disease_Type, 
           Disease_occurrence, 
           Sites_with_genus, 
           FOC), 
  file = "Diseases_FOC.csv"
)

# Print Diseases Prevalence by Site
write_csv(
  Diseases_prevalence_bySite %>% 
    select(Site, 
           Genus, 
           Disease_Type, 
           Count, 
           Transect_Area,
           Site_Count,
           SiteArea,
           Prevalence), 
  file = "Diseases_Prevalence_BySite.csv"
)

# Print Overall Diseases Prevalence
write_csv(
  Diseases_prevalence_overall %>%
    select(Genus, 
           Disease_Type, 
           Overall_Count_Disease, 
           Total_Disease_Area,
           Total_byGenera,
           Total_CoralCover_Area,
           Prevalence), 
  file = "Diseases_Prevalence_Overall.csv"
)

# Print Overall Genus Proportions
write_csv(CoralCover_GenusProportion, file = "Genus_Proportion_Overall.csv")

# Print coral densities and relative abundances by site
write_csv(CoralCover_RelAb, file = "Genus_Proportion_BySite.csv")

# Print substrate composition percentages
write_csv(PIT_summary %>%
            filter(Class_Percentage > 0),
          file = "Substrate_Composition.csv")

# Print size class proportions by site
write_csv(CoralCover_SizeClass_Proportions, file = "Proportion of size classes.csv")

# Print site-level disease prevalence
write_csv(Diseases_prevalence_SitesOnly, file = "Site-level disease prevalence.csv")

# Print site-level disease and coral counts (supplementary table)
Diseases_SiteCounts_Supplementary <- Diseases_prevalence_SitesOnly %>%
  mutate(Area_Ratio = Transect_Area / SiteArea) %>%
  select(Site, TotalDisease_Site, TotalCoralCount, Area_Ratio, Prevalence)

write_csv(Diseases_SiteCounts_Supplementary, file = "Site-level disease and coral counts.csv")

# Print coral density by genus and site (supplementary table, wide format)
CoralCover_GenusDensity_BySite <- CoralCover_filtered %>%
  group_by(Site, Genus_collapsed, SiteArea) %>%
  summarize(TotalCount = sum(TotalGenus), .groups = "drop") %>%
  mutate(Density = TotalCount / SiteArea,
         Site_Label = paste0(Site, " (", SiteArea, " m²)"))

CoralCover_Density_Table <- CoralCover_GenusDensity_BySite %>%
  select(Genus_collapsed, Site_Label, Density) %>%
  pivot_wider(names_from = Site_Label, values_from = Density, values_fill = 0) %>%
  arrange(Genus_collapsed)

write_csv(CoralCover_Density_Table, file = "Coral_Density_byGenus_and_Site.csv")
