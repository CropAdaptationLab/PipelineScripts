### ------------------------------------------------------------
### PCIL SELECTION PIPELINE — EXAMPLE WORKFLOW
### ------------------------------------------------------------

# Load functions from gists
source("https://gist.githubusercontent.com/claracruet/b6ade06ffa38c1e6bb97c813621632ea/raw/27fa2df7d1ed23ccc6faae951ac95e6426520fcf/load_pcil_data.R")
source("https://gist.githubusercontent.com/claracruet/189e3a4a2aabf0527ef0845832597439/raw/e360eaa231b57a31077358d2766dc423dac93b3a/select_pcil_positive.R")
source("https://gist.githubusercontent.com/claracruet/3f758a2f7d74a7d2f8278309b9500f67/raw/51d1952f4936e0f55a817474b3d8da8e83ce2e56/select_pcil_negative.R")


# ------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------

pcil_data <- load_pcil_data()

# Seed availability (defines "available world")
seed_metadata <- read.csv("~/Google Drive/My Drive/Clara lab folder/Program_Management/PCIL/PCIL_seed_source/CLARA_USE_ONLY_PCIL_MasterList_Seed_Harvested_2026-04-27.csv")

available_ids <- unique(seed_metadata$GENOTYPE.SAMPLE.ID)


# ------------------------------------------------------------
# INPUT SCENARIOS
# ------------------------------------------------------------

# Single gene
input_gene_single <- data.frame(Region = c("Sobic.002G215300"))

# Multiple genes (vector input)
input_gene_multi <- c("Sobic.002G215300", "Sobic.007G194300", "Sobic.010G205500")

# Explicit genomic regions
input_region <- data.frame(
  Region = c("Sobic.002G215300", "Sobic.007G194300", "Sobic.010G205500"),
  Chr    = c("Chr02", "Chr07", "Chr10"),
  Start  = c(62016113, 66054120, 56418423),
  End    = c(62018916, 66057692, 56423756)
)

# SNP / position input
input_snp_multi <- data.frame(
  Region = c("SNP1", "SNP2", "SNP3"),
  Chr    = c("Chr02", "Chr07", "Chr10"),
  pos    = c(1000, 2000, 3000)
)

input_snp_single <- data.frame(
  Region = "SNP1",
  Chr    = "Chr02",
  pos    = 1000
)


# ------------------------------------------------------------
# PCIL (+) SELECTION
# ------------------------------------------------------------

### Scenario 1 — Single gene

# Full population (genomics world)
pcil_pos_gene_full <- select_pcil_positive(
  pcil_data = pcil_data,
  input = input_gene_single,
  type = "gene"
)

# Available population (seed-constrained world)
pcil_pos_gene_available <- select_pcil_positive(
  pcil_data = pcil_data,
  input = input_gene_single,
  type = "gene",
  available_ids = available_ids
)


### Scenario 2 — Regions with selection

# Full population + top 2 per region
pcil_pos_region_full <- select_pcil_positive(
  pcil_data = pcil_data,
  input = input_region,
  type = "region",
  sel = 2
)

# Available population + top 2
pcil_pos_region_available <- select_pcil_positive(
  pcil_data = pcil_data,
  input = input_region,
  type = "region",
  sel = 2,
  available_ids = available_ids
)


# ------------------------------------------------------------
# PCIL (–) SELECTION
# ------------------------------------------------------------

### Scenario 1 — All PCIL (+)

# Best PCIL (–), full population
pcil_neg_gene_full <- select_pcil_negative(
  pcil_data = pcil_data,
  introgressions_in_region = pcil_pos_gene_full$pcil_positive,
  regions = pcil_pos_gene_full$regions
)

# Best PCIL (–), available population
pcil_neg_gene_available <- select_pcil_negative(
  pcil_data = pcil_data,
  introgressions_in_region = pcil_pos_gene_available$pcil_positive,
  regions = pcil_pos_gene_available$regions,
  available_ids = available_ids
)

# Top 4 PCIL (–) per PCIL (+)
pcil_neg_gene_top4 <- select_pcil_negative(
  pcil_data = pcil_data,
  introgressions_in_region = pcil_pos_gene_available$pcil_positive,
  regions = pcil_pos_gene_available$regions,
  n_neg = 4,
  available_ids = available_ids
)


### Scenario 2 — Using BEST PCIL (+)

# Best PCIL (–) for selected PCIL (+)
pcil_neg_bestpos_full <- select_pcil_negative(
  pcil_data = pcil_data,
  introgressions_in_region = pcil_pos_region_full$pcil_positive,
  pcil_positive_df = pcil_pos_region_full$best_lines,
  regions = pcil_pos_region_full$regions
)

# Available version
pcil_neg_bestpos_available <- select_pcil_negative(
  pcil_data = pcil_data,
  introgressions_in_region = pcil_pos_region_available$pcil_positive,
  pcil_positive_df = pcil_pos_region_available$best_lines,
  regions = pcil_pos_region_available$regions,
  available_ids = available_ids
)

# Top 4 PCIL (–) for best PCIL (+)
pcil_neg_bestpos_top4 <- select_pcil_negative(
  pcil_data = pcil_data,
  introgressions_in_region = pcil_pos_region_available$pcil_positive,
  pcil_positive_df = pcil_pos_region_available$best_lines,
  regions = pcil_pos_region_available$regions,
  n_neg = 4,
  available_ids = available_ids
)
