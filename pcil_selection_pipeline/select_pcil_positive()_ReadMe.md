### ------------------------------------------------------------
### 2) select_pcil_positive()
###   ARGUMENTS
###   ------------------------------------------------------------
###   pcil_data (list)
###     - Output from load_pcil_data()
###     - Contains all required datasets
###
###   type (character)
###     - One of: "gene", "region", "position"
###     - Controls how input is interpreted
###
###   input (varies)
###     - Defines genomic targets
###     - Depends on 'type':
###         A) type = "gene"
###           - character vector OR data.frame with Region column
###
###         B) type = "region"
###           - data.frame with:
###               Region | Chr | Start | End
###
###         C) type = "position"
###           - data.frame with:
###               Region | Chr | pos
###
###   sel (integer, optional)
###     - Number of top PCIL (+) lines to retain per region
###     - If NULL - all PCIL (+) are returned
###
###   donor_thresh (numeric, default = 0.75)
###     - Minimum mean donor fraction required
###
###   block_quantile (numeric, default = 0.75)
###     - Removes top fraction of largest introgression blocks
###
###   F_quantile (numeric, default = 0.25)
###     - Removes lowest fraction of inbreeding values
###
###   window (integer, default = 100)
###     - Used only for type = "position"
###     - Expands SNP/INDEL positions into regions
###
###   available_ids (character vector, optional)
###     - Vector of SampleIDs with available seed
###     - Restricts population BEFORE selection
###
###
###   ------------------------------------------------------------
###   KEY BEHAVIOR
###   ------------------------------------------------------------
###   - sel controls how many PCIL (+) are retained per region
###   - available_ids reduces population to those present in a list
###      *** such as those with available seed
###   - best_lines is ALWAYS a subset of pcil_positive
###
###   ------------------------------------------------------------
###   CORE STEPS
###   ------------------------------------------------------------
###     1. Convert input into standardized genomic regions
###     2. Identify introgression blocks fully covering each region
###     3. Merge genome-wide introgression metrics and inbreeding (F)
###
###   ------------------------------------------------------------
###   OPTIONAL SELECTION (sel = N)
###   ------------------------------------------------------------
###   If `sel` is provided, the function performs a structured
###   filtering and ranking process to select the top N PCIL (+)
###   lines per region.
###
###   Selection is applied independently per region.
###
###   Step 1 — Donor filter:
###     - Keep lines with sufficient donor signal
###       (mean_donor_frac ≥ donor_thresh)
###
###   Step 2 — Introgression size filter:
###     - Remove lines with largest introgression blocks
###       (top quantile defined by block_quantile)
###
###   Step 3 — Inbreeding filter:
###     - Remove least inbred lines
###       (bottom quantile defined by F_quantile)
###
###   Step 4 — Ranking:
###     Remaining lines ranked by:
###       1) lowest total introgressed Mb
###       2) fewest introgression blocks
###       3) highest inbreeding coefficient (F)
###
###   Step 5 — Selection:
###     - Retain top N lines per region (sel = N)
###
###   ------------------------------------------------------------
###   OUTPUT BEHAVIOR
###   ------------------------------------------------------------
###   If sel = NULL:
###     - Returns ALL PCIL (+) carriers
###
###   If sel is provided:
###     - Returns:
###         pcil_positive - ALL carriers
###         best_lines    - top N selected lines per region
###
###
###   ------------------------------------------------------------
###   AVAILABILITY MODE
###   ------------------------------------------------------------
###
###   available_ids:
###     - vector of SampleIDs
###     - restricts population BEFORE selection
###     - defines "operational" subset (lines with seed)
