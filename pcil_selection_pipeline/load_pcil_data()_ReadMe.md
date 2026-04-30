### 1) load_pcil_data()
###   ARGUMENTS
###   ------------------------------------------------------------
###
###   (No required arguments)
###
###   This function is designed to load pre-defined PCIL datasets
###   from external sources (e.g., GitHub gists or local files).
###
###   Output:
###     pcil_data - named list containing all required datasets
###
###   NOTE:
###     Paths or URLs may need to be edited depending on environment
###
###   OPERATION OF THE FUNCION
###   ------------------------------------------------------------
###   Loads and structures all required datasets:
###     - introgressions (per-line introgression blocks)
###     - genomewide_introgressions (total Mb, block counts)
###     - inbreeding_coefficient (F)
###     - gene_regions (gene coordinates)
###     - metadata (family, IDs)
###     - IBS_dis (pairwise genetic distance; 1 - IBS)
###
###   Output:
###     pcil_data (list containing all datasets)
