## Lolliplot Visualization of Variants in a Single Gene

The script [lolliplot_single_gene_command_line.R](https://github.com/CropAdaptationLab/Lab-Notebooks/blob/main/members/Tara/2026-projects/lolliplot_single_gene_command_line/lolliplot_single_gene_command_line.R) generates lollipop plots of SNPs and INDELs across a gene using the R package trackViewer. Variants are read from a VCF file output from SnpEff and plotted relative to gene features (gene body, exons, and UTRs) extracted from a [GFF3 annotation file](https://github.com/CropAdaptationLab/Lab-Notebooks/blob/main/members/Tara/2026-projects/lolliplot_single_gene_command_line/Sbicolorv5.1.gene_exons.gff3.gz).

### The script supports plotting:

* All variants in a gene

* List of selected variants in a gene

* A single variant in a gene


### The script requires the following R packages:

`trackViewer`

`VariantAnnotation`

`GenomicRanges`

`grid`

`gridExtra`

Install them in R if needed:

```
install.packages("BiocManager")
BiocManager::install(c("trackViewer","VariantAnnotation","GenomicRanges"))
install.packages("gridExtra")
```

### Required Files
1. Variant file (VCF)

A VCF containing SNPs and INDELs with SnpEff annotations in the ANN field. This file is provided as a command line input.

Example: [Sobic.003G241300_SNP_and_INDEL_merge.vcf](https://github.com/CropAdaptationLab/Lab-Notebooks/blob/main/members/Tara/2026-projects/lolliplot_single_gene_command_line/Sobic.003G241300_SNP_and_INDEL_merge.vcf)

The script extracts: chromosome, position, reference allele, alternate allele, and annotation

2. The compressed GFF3 file containing gene features for Sbicolor v5.1

File name: [Sbicolorv5.1.gene_exons.gff3.gz](https://github.com/CropAdaptationLab/Lab-Notebooks/blob/main/members/Tara/2026-projects/lolliplot_single_gene_command_line/Sbicolorv5.1.gene_exons.gff3.gz)

This GFF3 file needs to be in the same directory as working directory provided via the command line. You do not provide this file as input via command line.

### Script Usage

##### Option 1: To run the script from the command line:

`Rscript lolliplot_single_gene_command_line.R <vcf_file> <target_gene> <workdir> [variant_input]`

##### Option 2: To run the script from within RStudio:

In the console (bottom left pane of RStudio), set your variables and source the script:

```
commandArgs <- function(trailingOnly = TRUE) {
  c(
    "Sobic.003G241300_SNP_and_INDEL_merge.vcf",
    "Sobic.003G241300",
    "."
  )
}

source("lolliplot_single_gene_command_line.R")
```

#### Arguments:
- vcf_file:	Input VCF file
- target_gene:	Gene ID to visualize
- workdir:	Working directory containing input files
- variant:	Variant ID (required if not plotting all variants) in the format Chromosome:Position (ex: Chr03:64379137)
    - if you provide a list of variants, ensure the list has **NO** commas ex: `Chr03:64379137,Chr03:64379200,Chr03:64380004`

#### Example Commands
Plot all variants in a gene:

```
Rscript lolliplot_single_gene_command_line.R \
Sobic.003G241300_SNP_and_INDEL_merge.vcf \
Sobic.003G241300 \
.
```
**NOTE** : The default is to plot *ALL* variants in the gene

Output: [lolliplot_Sobic.003G241300_AllVariants.pdf](https://github.com/CropAdaptationLab/Lab-Notebooks/blob/main/members/Tara/2026-projects/lolliplot_single_gene_command_line/lolliplot_Sobic.003G241300_AllVariants.pdf)

This plot shows all SNPs and INDELs within the gene.

Plot a single variant in a gene: 

```
Rscript lolliplot_single_gene_command_line.R \
Sobic.003G241300_SNP_and_INDEL_merge.vcf \
Sobic.003G241300 \
. \
Chr03:64379137
```

Output: [lolliplot_Sobic.003G241300_SingleVariant.pdf](https://github.com/CropAdaptationLab/Lab-Notebooks/blob/main/members/Tara/2026-projects/lolliplot_single_gene_command_line/lolliplot_Sobic.003G241300_SingleVariants.pdf)

This plot shows only the specified variant within the gene.

Plot a list of variants in a gene: 

```
Rscript lolliplot_single_gene_command_line.R \
Sobic.003G241300_SNP_and_INDEL_merge.vcf \
Sobic.003G241300 \
. \
Chr03:64379137,Chr03:64379200,Chr03:64380004
```

Output: [lolliplot_Sobic.003G241300_SelectedVariants.pdf](https://github.com/CropAdaptationLab/Lab-Notebooks/blob/main/members/Tara/2026-projects/lolliplot_single_gene_command_line/lolliplot_Sobic.003G241300_SelectedVariants.pdf)

This plot shows only the specified variants within the gene.

### What is shown on the lollipop plot

- gene structure (gene body, exons, UTRs)

- variant positions

- variant annotations

- legend with variant annotation categories

- legend with gene feature types

- legend with SNP vs INDEL marker shapes
