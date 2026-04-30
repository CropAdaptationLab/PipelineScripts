### ------------------------------------------------------------
### 3) select_pcil_negative()
###
###   Identifies PCIL (–) lines (non-carriers) matched to PCIL (+)
###   lines for near-isogenic comparisons.
###
### ------------------------------------------------------------
###   ARGUMENTS
###   ------------------------------------------------------------
###   pcil_data (list)
###     - Output from load_pcil_data()
###
###   introgressions_in_region (data.frame)
###     - Output: pcil_positive from select_pcil_positive()
###     - Defines ALL PCIL (+) carriers
###
###   regions (data.frame)
###     - Output: regions from select_pcil_positive()
###     - Defines genomic intervals being evaluated
###
###   pcil_positive_df (data.frame, optional)
###     - Subset of PCIL (+) to process
###     - Typically: best_lines
###
###     - If NULL:
###         → all PCIL (+) are used
###
###     - If provided:
###         → only those PCIL (+) are matched
###
###   n_neg (integer, optional)
###     - Number of PCIL (–) matches per PCIL (+)
###
###     Behavior:
###       n_neg = NULL or 1:
###         → only best match returned
###
###       n_neg > 1:
###         → top N matches returned (ranked)
###
###   available_ids (character vector, optional)
###     - Restricts population BEFORE matching
###     - Ensures only seed-available lines are used
###
###
###   ------------------------------------------------------------
###   KEY BEHAVIOR
###   ------------------------------------------------------------
###
###   - Matching is always relative to PCIL (+)
###   - Same-family candidates are prioritized (soft constraint)
###   - IBS similarity is applied before genome refinement
###   - Some PCIL (+) may not yield valid matches and are skipped### 3) select_pcil_negative()
###  
###   Identifies PCIL (–) lines (non-carriers) matched to PCIL (+)
###   lines for near-isogenic comparisons.
###
###   ------------------------------------------------------------
###   CORE LOGIC
###   ------------------------------------------------------------
###
###   For each region and each PCIL (+):
###
###     STEP 1 — Candidate definition:
###       - All lines NOT carrying the region
###
###     STEP 2 — Family constraint (SOFT FILTER):
###       - Same-family candidates are prioritized
###       - If at least one same-family candidate exists:
###            - restrict selection to that family
###       - If no same-family candidates exist:
###            - fallback to full population
###
###       NOTE:
###         Same-family matching is a priority, NOT a requirement
###
###     STEP 3 — Tiered selection:
###
###       Tier 1 — IBS similarity:
###         - Compute IBS distance (1 - IBS)
###         - Select closest lines (genetically similar)
###
###       Tier 2 — Genome refinement:
###         - lowest total introgressed Mb
###         - fewest introgression blocks
###         - highest inbreeding coefficient (F)
###
###   ------------------------------------------------------------
###   SELECTION DETAILS
###   ------------------------------------------------------------
###
###   Tier 1 (IBS):
###     - If n_neg = 1:
###         - top 3 IBS candidates evaluated
###     - If n_neg > 1:
###         - top N IBS candidates retained
###
###   Tier 2 (refinement):
###     - Apply genome quality ranking
###     - Select final best and/or top N
###
###
###   ------------------------------------------------------------
###   AVAILABILITY MODE
###   ------------------------------------------------------------
###
###   available_ids:
###     - restricts population BEFORE matching
###     - ensures only deployable lines are considered
###
###
###   ------------------------------------------------------------
###   IMPORTANT NOTES
###   ------------------------------------------------------------
###
###   - Selection is performed independently per PCIL (+)
###   - Some PCIL (+) may not yield valid PCIL (–) due to:
###       * lack of candidates
###       * missing IBS data
###       * availability constraints
###
###     - these are skipped
#   These modes should be run separately and compared.
###
### ============================================================
