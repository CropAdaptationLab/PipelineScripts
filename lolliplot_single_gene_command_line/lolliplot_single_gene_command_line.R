library(trackViewer)
library(VariantAnnotation)
library(grid)
library(gridExtra)
library(GenomicRanges)

# ---------------------------
# Command line arguments
# ---------------------------
args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 3){
  stop("Usage: Rscript lolliplot.R <vcf_file> <target_gene> <workdir> [variant_list]")
}

vcf_file    <- args[1]
target_gene <- args[2]
workdir     <- args[3]

# Optional variant argument
variant <- if(length(args) >= 4) args[4] else NULL

# Set working directory
setwd(workdir)

# ---------------------------
# Load VCF
# ---------------------------
v1 <- readVcf(vcf_file, genome = "Sbicolor")

# ---------------------------
# Determine variants to plot
# ---------------------------
if(is.null(variant)){
  # Plot ALL variants
  vcf_ids <- rownames(v1)
  variant_ids_to_plot <- gsub("^(SNP|INDEL)_", "", vcf_ids)
  variant_ids_to_plot <- gsub("_", ":", variant_ids_to_plot, fixed=TRUE)
  
} else if(file.exists(variant)){
  # Read from file (one variant per line)
  variant_ids_to_plot <- readLines(variant)
  
} else {
  # Comma-separated list
  variant_ids_to_plot <- unlist(strsplit(variant, ","))
}

# ---------------------------
# Toggle color scheme
# ---------------------------
color_by <- "annotation"  # or "impact"

# ---------------------------
# Prepare GFF
# ---------------------------
gff.full <- read.table("Sbicolorv5.1.gene_exons.gff3.gz")
gff <- gff.full[gff.full$V3 %in% c("gene","exon","five_prime_UTR","three_prime_UTR"),]
gff$ID <- sub("^[^=]*=([^.]*\\.[^.]*).*", "\\1", gff$V9)
gff$len <- gff$V5 - gff$V4
gff$col <- ifelse(
  gff$V3=="exon","#ffdab9",
  ifelse(gff$V3 %in% c("five_prime_UTR","three_prime_UTR"),"grey90","grey40")
)
gff$height <- ifelse(
  gff$V3=="gene",0.01,
  ifelse(gff$V3 %in% c("exon","five_prime_UTR","three_prime_UTR"),0.07,0.03)
)

# ---------------------------
# Color maps
# ---------------------------
color_map <- c(
  "HIGH"="red","MODERATE"="orange",
  "LOW"="yellow","MODIFIER"="grey90"
)

annotation_color_map <- c(
  "splice_region_variant&intron_variant"="purple",
  "splice_region_variant"="purple",
  "intron_variant"="purple",
  "disruptive_inframe_insertion"="orangered",
  "frameshift_variant"="orangered",
  "frameshift_variant&stop_lost&splice_region_variant"="orangered",
  "frameshift_variant&stop_gained&splice_region_variant"="orangered",
  "frameshift_variant&stop_lost"="orangered",
  "splice_acceptor_variant&intron_variant"="orangered",
  "splice_donor_variant&intron_variant"="orangered",
  "splice_acceptor_variant&splice_region_variant&intron_variant"="orangered",
  "disruptive_inframe_deletion"="red4",
  "conservative_inframe_deletion"="red4",
  "frameshift_variant&start_lost"="red4",
  "missense_variant"="gold",
  "synonymous_variant"="yellowgreen",
  "splice_region_variant&synonymous_variant"="yellowgreen",
  "5_prime_UTR_variant"="lightblue",
  "5_prime_UTR_premature_start_codon_gain_variant"="darkblue",
  "3_prime_UTR_variant"="lightblue",
  "intron_variant"="grey90"
)

# ---------------------------
# Extract variant info
# ---------------------------
ann <- info(v1)$ANN
all_annotation <- sapply(ann, function(x) strsplit(as.character(x[1]), "\\|")[[1]][2])
all_impact     <- sapply(ann, function(x) strsplit(as.character(x[1]), "\\|")[[1]][3])
all_gene       <- sapply(ann, function(x) strsplit(as.character(x[1]), "\\|")[[1]][4])

all_ref <- as.character(ref(v1))
all_alt <- sapply(alt(v1), function(x) as.character(x[1]))
all_af <- (unlist(info(v1)$AF))*100
all_af <- ifelse(all_af > 50, 100 - all_af, all_af)

# ---------------------------
# Start from all variants
# ---------------------------
variants_gr <- rowRanges(v1)

# Keep only variants annotated to target gene
gene_idx <- all_gene == target_gene
variants_gr <- variants_gr[gene_idx]
all_annotation <- all_annotation[gene_idx]
all_impact <- all_impact[gene_idx]
all_ref <- all_ref[gene_idx]
all_alt <- all_alt[gene_idx]
all_af <- all_af[gene_idx]

# Variant IDs for filtering
variant_pos_id <- paste0(as.character(seqnames(variants_gr)), ":", start(variants_gr))

# Filter based on user input
keep_idx <- variant_pos_id %in% variant_ids_to_plot

# Only check for missing variants if user provided input
if(!is.null(variant)){
  missing_variants <- setdiff(variant_ids_to_plot, variant_pos_id)
  if(length(missing_variants) > 0){
    warning("These variants were not found in the gene/VCF: ",
            paste(missing_variants, collapse=", "))
  }
}

variants_gr <- variants_gr[keep_idx]
all_annotation <- all_annotation[keep_idx]
all_impact <- all_impact[keep_idx]
all_ref <- all_ref[keep_idx]
all_alt <- all_alt[keep_idx]
all_af <- all_af[keep_idx]

if(length(variants_gr)==0){
  stop("No matching variants found for the target gene.")
}

# Assign labels
names(variants_gr) <- paste0(variant_pos_id[keep_idx], " | [", all_ref, "/", all_alt, "] | ", all_annotation)

# Lollipop metadata
variants_gr$score <- all_af
variants_gr$color <- if(color_by=="impact") color_map[all_impact] else annotation_color_map[all_annotation]
variants_gr$SNPsideID <- "top"

# ---------------------------
# Prepare gene features
# ---------------------------
gff_gene <- gff[gff$ID == target_gene,]
chr <- unique(gff_gene$V1)[1]
gene_strand <- unique(gff_gene$V7[gff_gene$V3=="gene"])[1]

plot_range <- GRanges(chr, IRanges(min(gff_gene$V4)-2000, max(gff_gene$V5)+2000))

features_gr <- GRanges(chr, IRanges(gff_gene$V4, width=gff_gene$len),
                       strand = gene_strand,
                       fill = gff_gene$col,
                       height = gff_gene$height)

# ---------------------------
# Shape assignments
# ---------------------------
variants_gr$shape <- ifelse(width(variants_gr) == 1, "circle", "triangle_point_down")
end(variants_gr) <- start(variants_gr)

# ---------------------------
# Output file
# ---------------------------
if(is.null(variant)){
  file_type <- "AllVariants"
} else if(length(variant_ids_to_plot) == 1){
  file_type <- "SingleVariants"
} else {
  file_type <- "SelectedVariants"
}

pdf_file <- paste0("lolliplot_", target_gene, "_", file_type, ".pdf")

# ---------------------------
# Plot
# ---------------------------
pdf(pdf_file, width=14, height=12)
grid.newpage()

pushViewport(viewport(layout=grid.layout(nrow=2, ncol=2,
                                         widths=unit(c(0.82,0.18),"npc"),
                                         heights=unit(c(0.85,0.15),"npc"))))

pushViewport(viewport(layout.pos.row=1, layout.pos.col=1))

variants_gr_names <- paste0(as.character(seqnames(variants_gr)), ":", start(variants_gr))

vcf_variant_names <- rownames(v1)
vcf_variant_names <- gsub("^(SNP|INDEL)_", "", vcf_variant_names)
vcf_variant_names <- gsub("_", ":", vcf_variant_names, fixed=TRUE)

variant_idx <- match(variants_gr_names, vcf_variant_names)

gt <- geno(v1)$GT[variant_idx, , drop=FALSE]
gt_num <- matrix(as.numeric(gt), nrow=nrow(gt))
af <- (rowMeans(gt_num, na.rm=TRUE) / 2) * 100

mcols(variants_gr)$score <- af

lolliplot(variants_gr, features_gr, ranges=plot_range, jitter="node",
          ylab=paste(target_gene, "\nMAF (%)"),
          yaxis=c(0,25,50,75,100),
          legend=NULL, cex=0.8, srt=45, newpage=FALSE)

popViewport()

# ---------------------------
# Gene feature legend (bottom-left)
# ---------------------------
pushViewport(viewport(layout.pos.row=2, layout.pos.col=1))

feature_labels <- c("gene","exon","5'UTR","3'UTR")
feature_colors <- c("grey40","#ffdab9","grey90","grey90")

# Spread horizontally across the panel
x_positions <- seq(0.4, 0.6, length.out=length(feature_labels))
y_box <- 0.6
y_label <- 0.3
box_width <- 0.05
box_height <- 0.10

for(i in seq_along(feature_labels)){
  grid.rect(x=x_positions[i], y=y_box, width=box_width, height=box_height,
            gp=gpar(fill=feature_colors[i], col="black"))
  grid.text(feature_labels[i], x=x_positions[i], y=y_label,
            just="center", gp=gpar(fontsize=9))
}

popViewport()  # exit bottom-left panel

# ---------------------------
# Legends (unchanged)
# ---------------------------
pushViewport(viewport(layout.pos.row=1, layout.pos.col=2))

# Only annotations present in the plotted variants
unique_annotations <- unique(all_annotation)
legend_labels <- intersect(unique_annotations, names(annotation_color_map))
legend_colors <- annotation_color_map[legend_labels]

# Wrap long labels
wrapped_labels <- gsub("&", "&\n", legend_labels)
wrapped_labels <- gsub("premature_", "premature_\n", wrapped_labels)

n_labels <- length(legend_labels)

# Dynamically adjust the viewport height based on number of labels
legend_height <- min(0.80, 0.05 + 0.03*n_labels)  # max height 0.85 npc
legend_y_center <- 0.235 # keep centered vertically

# Compute y positions evenly spaced within the legend viewport
if(n_labels == 1){
  y_positions <- legend_y_center
} else {
  y_positions <- seq(from=legend_y_center + legend_height/2,
                     to=legend_y_center - legend_height/2,
                     length.out=n_labels)
}

# Draw legend title slightly above the top of legend
grid.text("Annotation", x=0.001, y=max(y_positions) + 0.05, just="left",
          gp=gpar(fontsize=10, fontface="bold"))

# Draw points and labels
for(i in seq_along(legend_labels)){
  grid.points(x=0.05, y=y_positions[i], pch=23,
              size=unit(0.035,"npc"),
              gp=gpar(col="black", fill=legend_colors[i]))
  grid.text(wrapped_labels[i], x=0.09, y=y_positions[i],
            just="left", gp=gpar(fontsize=8))
}

popViewport()

shape_legend <- legendGrob(labels=c("SNP","INDEL"),
                           pch=c(16,25),
                           gp=gpar(col="black", fontsize=8),
                           hgap=unit(0.4,"cm"),
                           vgap=unit(0.2,"cm"))

pushViewport(viewport(layout.pos.row=2, layout.pos.col=2))
pushViewport(viewport(x=0.15, y=1))
grid.draw(shape_legend)

popViewport(2)
popViewport()

dev.off()
