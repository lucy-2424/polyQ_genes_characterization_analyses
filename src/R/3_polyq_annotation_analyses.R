library(dplyr)
library(data.table) 
library(ghql)
library(ggplot2)
library(ggrepel)
library(biomaRt)
library(readxl)
library(janitor)
library(tidyverse)
library(viridisLite)
library(viridis)

# Load the polyq data 
l2g_combined <- fread("/research/2023_polyQ/data/polyQ_otg_29_09-23.tsv")
polyQ_disorder_genes <- read_excel("/research/2023_polyQ/data/polyQ_disorders_2023.xlsx",sheet = "Sheet1")

# Get ensemble polyq gene IDs
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", archive=FALSE, verbose=TRUE)
filters <- listFilters(ensembl)
attributes <- listAttributes(ensembl)
filterlist <- unique(polyQ_disorder_genes$Gene) 
polyq_gene_info <- getBM(attributes = c("hgnc_symbol","ensembl_gene_id", "start_position", "end_position"),filters = c("hgnc_symbol"),values = filterlist, mart = ensembl)

# Rename some column names
polyq_gene_info <- polyq_gene_info %>%
  rename(symbol = "hgnc_symbol",
         id = "ensembl_gene_id", 
         start_position = "start_position", 
         end_position = "end_position" )

# Add the start and end positions for each gene to the dataset with the polyq associations
l2g_combined <- left_join(l2g_combined, polyq_gene_info, by = c("id", "symbol"))

# Filter for the necessary columns that are required to retrieve the annotations 
l2g_combined_vep <- l2g_combined %>%
  dplyr::select(c(variant.chromosome, variant.position, variant.rsId, variant.refAllele,
                  variant.altAllele))

# Save the file to retrieve annotations from ensembl vep
write.table(l2g_combined_vep, "/research/2023_polyQ/data/l2g_combined_no_annot.txt", quote = F, row.names = F)

## Retrieve annotations using VEP in bash
# Load annotated data
l2g_combined_annotated <- fread("/research/2023_polyQ/data/l2g_combined_annot.txt")

# Select on the required columns 
l2g_combined_annotated <- l2g_combined_annotated %>%
  dplyr::select(c(`#Uploaded_variation`, Location, Consequence))

# Rename columns 
l2g_combined_annotated <- l2g_combined_annotated %>%
  rename(variant.rsId = `#Uploaded_variation`,
         location = Location,
         consequence = Consequence) 

# Retain unique rows
l2g_combined_annotated <- unique(l2g_combined_annotated)

# Create a new column called location in the l2g_combined
l2g_combined$location <- paste(l2g_combined$variant.chromosome, l2g_combined$variant.position, sep = ":")

# Combine the annotated data frame with the original L2G 
l2g_combined_merge <- merge(l2g_combined, l2g_combined_annotated, by = "location")

# Check for missing information
# Check if all the missing rows have count phenotype, l2g<0.5 or missing pmid
missing_rows <- anti_join(l2g_combined, l2g_combined_merge)
missing_rows <- missing_rows %>%
  filter(L2G>0.5)

missing_rows <- missing_rows %>% 
  filter(!(study.traitReported %like% " cell")) %>% 
  filter(!(study.traitReported %like% " count")) %>%
  filter(!(study.traitReported %like% " level")) %>%
  filter(!(study.traitReported %like% " fraction")) %>%
  filter(!(study.traitReported %like% " pressure")) %>%
  filter(!(study.traitReported %like% "Mean ")) %>%
  filter(!(study.traitReported %like% "Lung ")) %>%
  filter(!(study.traitReported %like% "Protein quant")) %>%
  filter(!(study.traitReported %like% "Height")) %>%
  filter(!(study.traitReported %like% " percent")) %>%
  filter(!(study.traitReported %like% "Waist")) %>%
  filter(!(study.traitReported %like% "Platelet")) %>%
  filter(!(study.traitReported %like% "crit")) %>%
  filter(!(study.traitReported %like% " density")) %>%
  filter(!(study.traitReported %like% "Hemoglobin")) %>%
  filter(!(study.traitReported %like% "filtration rate")) %>%
  filter(!(study.traitReported %like% "Spleen volume")) %>%
  filter(!is.na(study.pmid))

# One missing row passing the qc that needs to be added to the overall data set
missing_rows$consequence <- c("intron_variant")
l2g_combined_merge <- l2g_combined_merge %>% dplyr::select(!c(variant.rsId.y))
l2g_combined_merge <- l2g_combined_merge %>% rename(variant.rsId = variant.rsId.x)
l2g_combined_merge <- rbind(l2g_combined_merge, missing_rows)

# Re-categorize the consequence column 
l2g_combined_merge$consequence <- gsub("_", " ", l2g_combined_merge$consequence)
l2g_combined_merge$consequence <- gsub("5 prime UTR variant", "UTR variant", l2g_combined_merge$consequence)
l2g_combined_merge$consequence <- gsub("3 prime UTR variant", "UTR variant", l2g_combined_merge$consequence)
l2g_combined_merge$consequence <- gsub("splice acceptor variant", "splice variant", l2g_combined_merge$consequence)
l2g_combined_merge$consequence <- gsub("splice donor region variant", "splice variant", l2g_combined_merge$consequence)
l2g_combined_merge$consequence <- gsub("splice region variant", "splice variant", l2g_combined_merge$consequence)
l2g_combined_merge$consequence <- gsub("splice polypyrimidine tract variant", "splice variant", l2g_combined_merge$consequence)
l2g_combined_merge$consequence <- gsub("upstream gene variant", "other variant", l2g_combined_merge$consequence)
l2g_combined_merge$consequence <- gsub("downstream gene variant", "other variant", l2g_combined_merge$consequence)
l2g_combined_merge$consequence <- gsub("non coding transcript exon variant", "other variant", l2g_combined_merge$consequence)

# Generate data sets for just the significant genes for plotting later
htt <- l2g_combined_merge %>%
  dplyr::filter(symbol == "HTT")
atxn2 <- l2g_combined_merge %>%
  dplyr::filter(symbol == "ATXN2")
atxn1 <- l2g_combined_merge %>%
  dplyr::filter(symbol == "ATXN1")
atxn7 <- l2g_combined_merge %>%
  dplyr::filter(symbol == "ATXN7")
cacna1a <- l2g_combined_merge %>%
  dplyr::filter(symbol == "CACNA1A")

# Reformat and take out the excess wording from the traits
l2g_combined_merge$study.traitReported <- gsub(" \\[EA\\])", "", l2g_combined_merge$study.traitReported)
l2g_combined_merge$study.traitReported <- gsub(" \\(years of education\\)", "", l2g_combined_merge$study.traitReported)
l2g_combined_merge$study.traitReported <- gsub(" \\(special factor of neuroticism\\)", "", l2g_combined_merge$study.traitReported)
l2g_combined_merge$study.traitReported <- gsub(" \\(healthspan, parental lifespan or longevity)", "", l2g_combined_merge$study.traitReported)
l2g_combined_merge$study.traitReported <- gsub(" \\(multivariate analysis\\)", "", l2g_combined_merge$study.traitReported)
l2g_combined_merge$study.traitReported <- gsub(" \\[MTAG\\]", "", l2g_combined_merge$study.traitReported)
l2g_combined_merge$study.traitReported <- gsub(" \\(MTAG\\)", "", l2g_combined_merge$study.traitReported)
l2g_combined_merge$study.traitReported <- gsub(" \\(Hounsfield unit scale\\)", "", l2g_combined_merge$study.traitReported)
l2g_combined_merge$study.traitReported <- gsub(" \\(pleiotropy\\)", "", l2g_combined_merge$study.traitReported)
l2g_combined_merge$study.traitReported <- gsub(" \\(years of education\\)", "", l2g_combined_merge$study.traitReported)
l2g_combined_merge$study.traitReported <- gsub(" \\[EA\\]", "", l2g_combined_merge$study.traitReported)
l2g_combined_merge$study.traitReported <- gsub(" \\[conditional-joint\\]", "", l2g_combined_merge$study.traitReported)
l2g_combined_merge$study.traitReported <- gsub(" \\[MTAG\\]", "", l2g_combined_merge$study.traitReported)
l2g_combined_merge$study.traitReported <- gsub(" \\(MTAG\\)", "", l2g_combined_merge$study.traitReported)
l2g_combined_merge$study.traitReported <- gsub("/atrial flutter", "", l2g_combined_merge$study.traitReported)

## Create functions for the repetitive operations 
# Function to remove duplicate associations 
rank_and_deduplicate <- function(data) {
  data %>%
    group_by(study.traitReported) %>%
    slice_max(order_by = L2G, n = 1) %>%
    ungroup()
}

# Function to select relevant associations (remove count based, etc)
filter_phenotypes <- function(data) {
  data %>%
    filter(!(study.traitReported %like% " cell")) %>%
    filter(!(study.traitReported %like% " count")) %>%
    filter(!(study.traitReported %like% " level")) %>%
    filter(!(study.traitReported %like% " fraction")) %>%
    filter(!(study.traitReported %like% " pressure")) %>%
    filter(!(study.traitReported %like% "Mean ")) %>%
    filter(!(study.traitReported %like% "Lung ")) %>%
    filter(!(study.traitReported %like% "Protein quant")) %>%
    filter(!(study.traitReported %like% "Height")) %>%
    filter(!(study.traitReported %like% " percent")) %>%
    filter(!(study.traitReported %like% "Waist")) %>%
    filter(!(study.traitReported %like% "Platelet")) %>%
    filter(!(study.traitReported %like% "crit")) %>%
    filter(!(study.traitReported %like% " density")) %>%
    filter(!(study.traitReported %like% "Hemoglobin")) %>%
    filter(!(study.traitReported %like% "filtration rate")) %>%
    filter(!(study.traitReported %like% "Spleen volume")) %>%
    filter(!is.na(study.pmid)) %>%
    filter(L2G > 0.5)
}

# Function to plot significant associations within gene area shaded
plot_gene <- function(data, xmin, xmax, gene_name, chr_position, output_path, data_pheno, position) {
  ggplot(data, aes(y = L2G, x = variant.position, col = consequence)) +
    geom_rect(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, 
              color = NA, fill = 'grey80', alpha = 0.2) +
    geom_point(size = 1) + theme_bw() +
    theme(legend.position = "right", 
      axis.title = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(size = 10, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(size = 10), 
      axis.title.x = element_text(size = 8),  
      axis.title.y = element_text(size = 8)) +
    ylim(c(0, 1)) + xlab(chr_position) + ylab("L2G pipeline score") +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
    geom_text_repel(
      data = filter(data_pheno , L2G > 0.5), 
      aes(label = study.traitReported),
      max.overlaps = Inf, box.padding = 0.5,
      point.padding = 0.5, segment.color = "grey50", 
      size = 3.5) +
    scale_color_brewer(palette = "Dark2") + 
    annotate(
      "text", 
      x = position, y = 1.00, 
      label = bquote(italic(.(gene_name))), 
      size = 4, color = "black", 
      family = "serif"
    )
  
  ggsave(output_path, height = 4.36, width = 7.36, units = 'in')
}

# Apply the filtering and ranking functions for relevant genes
htt_pheno <- rank_and_deduplicate(filter_phenotypes(l2g_combined_merge %>% filter(symbol == "HTT")))
atxn1_pheno <- rank_and_deduplicate(filter_phenotypes(l2g_combined_merge %>% filter(symbol == "ATXN1")))
atxn2_pheno <- rank_and_deduplicate(filter_phenotypes(l2g_combined_merge %>% filter(symbol == "ATXN2")))
atxn7_pheno <- rank_and_deduplicate(filter_phenotypes(l2g_combined_merge %>% filter(symbol == "ATXN7")))
cacna1a_pheno <- rank_and_deduplicate(filter_phenotypes(l2g_combined_merge %>% filter(symbol == "CACNA1A")))

# Remove one duplicate association having the same L2G score in atxn1_pheno 
atxn1_pheno <- atxn1_pheno %>%
  group_by(study.traitReported) %>%
  slice_head(n = 1) %>%
  ungroup()

# Apply the plot function to generate plots for each gene
plot_gene(
  data = htt, 
  data_pheno = htt_pheno,
  xmin = 3041363, xmax = 3243957, 
  position = 3150000,
  gene_name = "HTT", 
  chr_position = "Chromosome 4 position (GRCh38)", 
  output_path = "/research/2023_polyQ/results/polyq_genes_htt_annotated_associations.pdf"
)

plot_gene(
  data = atxn1,
  data_pheno = atxn1_pheno,
  xmin = 16299112, xmax = 16761491, 
  position = 16400000,
  gene_name = "ATXN1", 
  chr_position = "Chromosome 6 position (GRCh38)", 
  output_path = "/research/2023_polyQ/results/polyq_genes_atxn1_annotated_associations.pdf"
)

plot_gene(
  data = atxn2, 
  data_pheno = atxn2_pheno,
  xmin = 111443485, xmax = 111599676, 
  position = 111500000,
  gene_name = "ATXN2", 
  chr_position = "Chromosome 12 position (GRCh38)", 
  output_path = "/research/2023_polyQ/results/polyq_genes_atxn2_annotated_associations.pdf"
)

plot_gene(
  data = atxn7, 
  data_pheno = atxn7_pheno,
  position = 63940000,
  xmin = 63863155, xmax = 64003462, 
  gene_name = "ATXN7", 
  chr_position = "Chromosome 3 position (GRCh38)", 
  output_path = "/research/2023_polyQ/results/polyq_genes_atxn7_annotated_associations.pdf"
)

plot_gene(
  data = cacna1a, 
  data_pheno = cacna1a_pheno,
  xmin = 13206442, xmax = 13624489, 
  position = 13400000,
  gene_name = "CACNA1A", 
  chr_position = "Chromosome 19 position (GRCh38)", 
  output_path = "/research/2023_polyQ/results/polyq_genes_cacna1a_annotated_associations.pdf"
)
