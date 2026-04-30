######## DEVELOPMENT/EDIT OF DARTSEQ TO VCF FILES #############

# I copied the vcf files from Tara to my directory
Tara's files:
vcf1-/pl/active/Morris_CSU/PCIL/updated_DArTseq_VCFs/OrderAppendix_1_DS25-10850/Report_DS25-10850_1_moreOrders_SNP_vcf_2.vcf.gz
vcf2-/pl/active/Morris_CSU/PCIL/updated_DArTseq_VCFs/OrderAppendix_1_DS25-10850/Report_DS25-10850_1_moreOrders_SNP_vcf_3.vcf.gz

# location im working from
vcf1- '/pl/active/Morris_CSU/Clara_Cruet/Report_DS25-10850_1_moreOrders_SNP_vcf_2.vcf.gz'
vcf2- '/pl/active/Morris_CSU/Clara_Cruet/Report_DS25-10850_1_moreOrders_SNP_vcf_3.vcf.gz'

#creating header
cat <<EOF > header.txt
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=Chr01>
##contig=<ID=Chr02>
##contig=<ID=Chr03>
##contig=<ID=Chr04>
##contig=<ID=Chr05>
##contig=<ID=Chr06>
##contig=<ID=Chr07>
##contig=<ID=Chr08>
##contig=<ID=Chr09>
##contig=<ID=Chr10>
##contig=<ID=UNK>
##contig=<ID=scaffold_934>
##contig=<ID=scaffold_2382>
EOF

#append samples
zcat Report_DS25-10850_1_moreOrders_SNP_vcf_2.vcf.gz | grep "^#CHROM" >> header.txt

bcftools reheader \
  -h header.txt \
  -o Report_DS25-10850_1_moreOrders_SNP_vcf_2_reheader.vcf.gz \
  Report_DS25-10850_1_moreOrders_SNP_vcf_2.vcf.gz  
bcftools sort Report_DS25-10850_1_moreOrders_SNP_vcf_2_reheader.vcf.gz -Oz -o Report_DS25-10850_1_moreOrders_SNP_vcf_2_reheader_sorted.vcf.gz
bcftools index Report_DS25-10850_1_moreOrders_SNP_vcf_2_reheader_sorted.vcf.gz


#creating header
cat <<EOF > header2.txt
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=Chr01>
##contig=<ID=Chr02>
##contig=<ID=Chr03>
##contig=<ID=Chr04>
##contig=<ID=Chr05>
##contig=<ID=Chr06>
##contig=<ID=Chr07>
##contig=<ID=Chr08>
##contig=<ID=Chr09>
##contig=<ID=Chr10>
##contig=<ID=UNK>
##contig=<ID=scaffold_265>
##contig=<ID=scaffold_934>
##contig=<ID=scaffold_1040>
##contig=<ID=scaffold_1040:24938>
##contig=<ID=scaffold_1469>
##contig=<ID=scaffold_1469:34793>
##contig=<ID=scaffold_2382>
##contig=<ID=scaffold_2379>
##contig=<ID=scaffold_2379:4128>
##contig=<ID=scaffold_2933>
##contig=<ID=scaffold_2933:16476>
EOF

#append samples
zcat Report_DS25-10850_1_moreOrders_SNP_vcf_3.vcf.gz | grep "^#CHROM" >> header2.txt

bcftools reheader \
  -h header2.txt \
  -o Report_DS25-10850_1_moreOrders_SNP_vcf_3_reheader.vcf.gz \
  Report_DS25-10850_1_moreOrders_SNP_vcf_3.vcf.gz  
bcftools sort Report_DS25-10850_1_moreOrders_SNP_vcf_3_reheader.vcf.gz -Oz -o Report_DS25-10850_1_moreOrders_SNP_vcf_3_reheader_sorted.vcf.gz
bcftools index Report_DS25-10850_1_moreOrders_SNP_vcf_3_reheader_sorted.vcf.gz


# bind the files
bcftools concat \
  -a \
  -Oz \
  -o Report_DS25-10850_1_moreOrders_SNP_vcf_merged_reheader_sorted.vcf.gz \
  Report_DS25-10850_1_moreOrders_SNP_vcf_2_reheader_sorted.vcf.gz \
  Report_DS25-10850_1_moreOrders_SNP_vcf_3_reheader_sorted.vcf.gz

# final vcf file for dartseq
Report_DS25-10850_1_moreOrders_SNP_vcf_merged_reheader_sorted.vcf.gz

######## HETEROZYGOCITY/INBREEDING CALCULATIONS #############

# I'm going to compute genome wide het calculation in order to select the best PCIL lines

# I'm starting with this file
Report_DS25-10850_1_moreOrders_SNP_vcf_merged_reheader_sorted.vcf.gz

# getting stats fisrt
/projects/claramcb@colostate.edu/software/plink --vcf Report_DS25-10850_1_moreOrders_SNP_vcf_merged_reheader_sorted.vcf.gz \
      --double-id \
      --het \
      --allow-extra-chr \
      --out Report_DS25-10850_1_moreOrders_SNP_vcf_merged_reheader_sorted_hetero.txt


# outcome
Report_DS25-10850_1_moreOrders_SNP_vcf_merged_reheader_sorted_hetero.txt.het

# I will use the inbreeding coeficient rather than het moving forward
<img width="1242" height="198" alt="image" src="https://github.com/user-attachments/assets/d3ee156b-69b9-4bd1-9bf9-dd3722a22011" />

This is F coeficcient

############### IBS distance matrix calculation ###################

# calculating IBS distance matrix 1-ibs
/projects/claramcb@colostate.edu/software/plink --vcf Report_DS25-10850_1_moreOrders_SNP_vcf_merged_reheader_sorted.vcf.gz \
      --double-id \
      --allow-extra-chr \
      --maf 0.05 \
      --geno 0.1 \
      --distance square 1-ibs \
      --out pcil_merged_dartseq_ibs_maf0.05_mis0.1


# outputs for the pipeline:
pcil_merged_dartseq_ibs_maf0.05_mis0.1.log
pcil_merged_dartseq_ibs_maf0.05_mis0.1.mdist
pcil_merged_dartseq_ibs_maf0.05_mis0.1.mdist.id
pcil_merged_dartseq_ibs_maf0.05_mis0.1.nosex
