# Load relevant libraries
library(dplyr)
library(data.table) 
library(ghql)
library(jsonlite) 
library(ggplot2)
library(ggrepel)
library(forcats)
library(readxl)
library(biomaRt)
library(tidyverse)
library(janitor)
library(rrapply)

# Create the gene dataset
polyq_genes <- c("ATN1", "AR", "ATXN1", "ATXN2", "ATXN3", "CACNA1A", "ATXN7", "TBP", "THAP11", "HTT")

# Read in gnomAD information
gnomad <- fread("/research/references/gnomad/gnomad.v2.1.1.lof_metrics.by_gene.txt")
gnomad_select <- gnomad %>% dplyr::select(gene, oe_lof, pLI) 
polyq_genes_info <- dplyr::filter(gnomad_select, gene %in% polyq_genes)

# Plot the gnomAD information
polyq_genes_info %>%
  ggplot(aes(x = -oe_lof, y = pLI)) +
  geom_point(aes(color = gene), size = 3) + theme_bw() +
  scale_color_manual(values=c("THAP11" = "#addc30","TBP" ="#fde725","ATN1" = "#5ec962", "ATXN7" = "#28ae80", "ATXN3" = "#21918c","ATXN1" = "#2c728e",
                             "CACNA1A" = "#3b528b","HTT" = "#472d7b","ATXN2" = "#440154"))+
  geom_hline(yintercept = 0.9, linetype = "dotted", color = "red") + 
  annotate("text", x = -0.3, y = 0.9, label = "pLI = 0.9", vjust = -1, color = "red") +  
  ggrepel::geom_text_repel(aes(label = gene, color = gene), size = 3, fontface = "bold.italic") + 
  labs(x = "gnomAD -(observed/expected ratio LoF)", y = "gnomADpLI") +
  theme(legend.position = "none",
        axis.text = element_text(size = 12, face = "bold")) 
ggsave("/research/2023_polyQ/results/polyq_genes_gnomad.pdf", height=4.50, width=5.36, units='in')

# Read in protein atlas information
protein_atlas <- fread("/research/references/protein_atlas/proteinatlas.tsv", na.strings = "")
protein_atlas <- protein_atlas %>% 
  clean_names()
protein_atlas_select <- protein_atlas %>%
  dplyr::select(gene, rna_tissue_specificity, rna_single_cell_type_specificity)
polyq_genes_info <- merge(polyq_genes_info, protein_atlas_select, by = "gene")

# Plot RNA tissue specificity information for the genes
# Count occurrences of each gene within each RNA tissue specificity group
polyq_genes_info_summary <- polyq_genes_info %>%
  dplyr::group_by(rna_tissue_specificity, gene) %>%
  dplyr::summarize(count = n(), .groups = 'drop')

# Create the stacked bar plot
polyq_genes_info_summary %>% 
  ggplot(aes(x = rna_tissue_specificity, y = count, fill = rna_tissue_specificity)) +
  geom_col(position = "stack", color = "white", width = 1) +
  geom_text(aes(label = gene), 
            position = position_stack(vjust = 0.5), 
            fontface = "italic", 
            color = "black") +
  scale_fill_manual(values = c( "#189AB4", "#75E6DA", "#D4F1F4")) +
  theme_minimal() +
  theme(axis.title.y = element_blank(),  # Remove y-axis title
        axis.text.y = element_blank(),  # Remove y-axis text
        axis.ticks.y = element_blank(),  # Remove y-axis ticks
        panel.grid = element_blank(), # Remove the gid
        legend.position = "none") +  # Remove legend
  labs(x = "RNA Tissue Specificity", y = "")
ggsave("/research/2023_polyQ/results/polyq_genes_tissue_specificity.pdf", height=5.58, width=6.58, units='in')
rm(polyq_genes_info_summary)

# Assess and add interaction data from Open Targets
# Subset high quality (mi score greater than 0.42, Threshold based on https://doi.org/10.1101/2023.02.07.23285407)
# Retrieve polyQ gene ensembl ids
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", archive=FALSE, verbose=TRUE)
filters <- listFilters(ensembl)
attributes <- listAttributes(ensembl)
filterlist <- unique(polyq_genes_info$gene) 
emsembl_gene_info <- getBM(attributes = c("hgnc_symbol","ensembl_gene_id"),filters = c("hgnc_symbol"),values = filterlist, mart = ensembl)
polyq_genes_mi <- merge(polyq_genes_info, emsembl_gene_info, by.x="gene", by.y="hgnc_symbol", all.x=T)

# Set up to query Open Targets Platform API
otp_cli <- GraphqlClient$new(url = 'https://api.platform.opentargets.org/api/v4/graphql')
otp_qry <- Query$new()
otp_qry$query('intact_query', 'query InteractionsSectionQuery(
  $ensgId: String!
  $sourceDatabase: String
  $index: Int = 0
  $size: Int = 10
) {
  target(ensemblId: $ensgId) {
    id
    approvedName
    approvedSymbol
    interactions(
      sourceDatabase: $sourceDatabase
      page: { index: $index, size: $size }
    ) {
      count
      rows {
        intA
        intABiologicalRole
        targetA {
          id
          approvedSymbol
        }
        speciesA {
          mnemonic
        }
        intB
        intBBiologicalRole
        targetB {
          id
          approvedSymbol
        }
        speciesB {
          mnemonic
        }
        score
        count
        sourceDatabase
        evidences {
          evidenceScore
          hostOrganismScientificName
          interactionDetectionMethodMiIdentifier
          interactionDetectionMethodShortName
          interactionIdentifier
          interactionTypeShortName
          participantDetectionMethodA {
            miIdentifier
            shortName
          }
          participantDetectionMethodB {
            miIdentifier
            shortName
          }
          interactionDetectionMethodShortName
          expansionMethodShortName
          pubmedId
        }
      }
    }
  }
}')

# Execute the query 
variables <- list(ensgId= c("ENSG00000066427, ENSG00000111676, ENSG00000112592,
                            ENSG00000124788, ENSG00000141837, ENSG00000163635,
                            ENSG00000168286, ENSG00000169083, ENSG00000197386,
                            ENSG00000204842"),
                  sourceDatabase= c("intact"),
                  index= c(Int=0),
                  size= c(Int=1000))
(otp_cli$exec(otp_qry$queries$intact_query, variables, flatten = TRUE))

# JSON arrays automatically convert list into more specific R class using the function "fromJSON"
result <- fromJSON(otp_cli$exec(otp_qry$queries$intact_query, variables, flatten = TRUE))$data
result_df <- as.data.frame(result) %>% flatten()

# Write a function that fetches intact information for genes
# Use this later to automate as you can only fetch one gene at a time
fetch_intact <- function(current_ensgId) {
  variables = list(ensgId= current_ensgId,
                   sourceDatabase= c("intact"),
                   index= c(Int=0),
                   size= c(Int=1000))
  result_gene <- fromJSON(otp_cli$exec(otp_qry$queries$intact_query, variables, flatten = TRUE))$data
  if (!is.null(result_gene[["target"]][["interactions"]])){
    #result_gene_df <- lapply(result_gene,as.data.frame)
    result_gene_df <- as.data.frame(result_gene) %>% dplyr::select(target.interactions.rows.score, target.approvedSymbol)
    return(result_gene_df)
  }
}

# Extract intact info for all genes of interest
query_genes <- polyq_genes_info$ensembl_gene_id
query_genes <- na.omit(query_genes)

intact_combined <- NULL
for (hgnc_symbol in query_genes) {
  print(hgnc_symbol)
  temp_intact <- fetch_intact(hgnc_symbol)
  print("binding")
  intact_combined <- rbind(intact_combined, temp_intact)
}


# Incase you want to download all information without intializing rows otherwise skip
for (hgnc_symbol in query_genes) {
  print(hgnc_symbol)
  temp_intact <- fetch_intact(hgnc_symbol)
  
  print("binding")
  
  # Use dplyr's bind_rows to avoid row name issues
  intact_combined <- bind_rows(intact_combined, temp_intact)
}
otg_interactors <- as.data.frame(intact_combined)

# Clean the data and retain the relevant information
otg_interactors <- otg_interactors %>%
  dplyr::select(c(target.interactions.rows.targetB, target.interactions.rows.score, 
                  target.interactions.rows.intABiologicalRole, target.interactions.rows.intBBiologicalRole,
                  target.interactions.rows.count, target.approvedSymbol))

otg_interactors$target.interactions.rows.intABiologicalRole <- paste0("A:", otg_interactors$target.interactions.rows.intABiologicalRole)
otg_interactors$target.interactions.rows.intBBiologicalRole <- paste0("B:", otg_interactors$target.interactions.rows.intBBiologicalRole)

otg_interactors$biologicalRole <- paste(otg_interactors$target.interactions.rows.intABiologicalRole,
                                        otg_interactors$target.interactions.rows.intBBiologicalRole,
                                        sep = ", ")
rename(data, new_name = old_name)
otg_interactors <- otg_interactors %>%
  rename(interactorB = target.interactions.rows.targetB,
         score = target.interactions.rows.score,
         interactionEvidenceEntries = target.interactions.rows.count,
         gene = target.approvedSymbol) %>%
  dplyr::select(!c(target.interactions.rows.intABiologicalRole, target.interactions.rows.intBBiologicalRole))

# Save the reference data 
# Data downloaded 2024-10-23
write.table(otg_interactors, file = "/research/2023_polyQ/otg/data/polyQ_genes_molecular_interactions_interactors.txt", sep = "\t", row.names = FALSE)

otg_interactors <- fread("/research/2023_polyQ/data/polyQ_genes_molecular_interactions_interactors.txt")

polyq_int_high_q_count <- otg_interactors %>%
  dplyr::filter(score>0.42) %>% 
  group_by(gene) %>%
  summarise(int_0.42_count = n())
polyq_int_count <- otg_interactors %>%
  group_by(gene) %>%
  summarise(int_all_count = n())

polyq_genes_info <- merge(polyq_genes_info, polyq_int_high_q_count, by="gene")
polyq_genes_info <- merge(polyq_genes_info, polyq_int_count, by="gene")
polyq_genes_info <- polyq_genes_info %>%
  mutate(int_0.42_count = tidyr::replace_na(int_0.42_count, 0)) %>%
  mutate(int_all_count = tidyr::replace_na(int_all_count, 0))


# Plot molecular interaction information
polyq_int_high_q_count %>%
  mutate(name = fct_reorder(gene, int_0.42_count)) %>%
  ggplot(aes(x=name, y=int_0.42_count, fill = gene)) +
  geom_bar(stat="identity", width=.9) +
  geom_text(aes(label=int_0.42_count),
            hjust = -0.1, 
            size = 3,
            fontface = "bold") +
  labs(x = "", y = "High quality interactions (Intact MI score > 0.42)") +
  theme_minimal() +
  geom_hline(yintercept = 10, linetype = "dotted", color = "black") +  # Add vertical line
  annotate("text", x = 11.9, y = 173, size = 3,
           label = "10 interactions", angle = 15, vjust = 1, hjust = 1, color = "black", fontface = "bold") +
  theme(axis.text.y = element_text(face = "bold.italic"),
        axis.text.x = element_text(face = "bold"),
        legend.position = "none") +
  scale_fill_manual(values = c(
    "THAP11" = "#addc30", "TBP" = "#fde725", "ATN1" = "#5ec962",
    "ATXN7" = "#28ae80", "ATXN3" = "#21918c", "ATXN1" = "#2c728e",
    "CACNA1A" = "#3b528b", "HTT" = "#472d7b", "ATXN2" = "#440154"
  )) +
  coord_flip() 
ggsave("/research/2023_polyQ/results/polyq_genes_interactors.pdf", height=4.36, width=6.36, units='in')

# Read in Druggable Genome Information (https://www.dgidb.org/, analyzed 23 July 2024)
dgib <- fread("/research/2023_polyQ/data/categories_dgidb_2024_07_23.tsv")

# Change the column names 
dgib <- dgib %>%
  rename(gene = name, category = `name-2`)

# Remove clinically actionable as this appears related to cancer
dgib <- dgib %>% 
  dplyr::filter(category != "CLINICALLY ACTIONABLE")
table(dgib$category)
dgib_drug_genes <- dgib %>% 
  dplyr::filter(category=="DRUGGABLE GENOME") %>%
  .$gene
dgib_other_drug_genes <- dgib %>% 
  dplyr::filter(category!="DRUGGABLE GENOME") %>%
  .$gene
dgib_other_drug_genes <- unique(dgib_other_drug_genes)

# Add information to main data frame
polyq_genes_info$druggable_genome <- ifelse(polyq_genes_info$gene %in% dgib_drug_genes, "Yes", "No")
polyq_genes_info$druggable_other <- ifelse(polyq_genes_info$gene %in% dgib_other_drug_genes, "Yes", "No")
polyq_genes_info$druggable_any <- ifelse(polyq_genes_info$druggable_genome=="Yes"| polyq_genes_info$druggable_other=="Yes", "Yes", "No")

# Base metrics on findings from https://www.nature.com/articles/s41588-024-01854-z
# https://www.nature.com/articles/s41588-024-01854-z
polyq_genes_info <- mutate(polyq_genes_info, pLI = ifelse(pLI < 0.9, "Tolerant", "Intolerant"))
polyq_genes_info$interacting_partners <- ifelse(polyq_genes_info$int_0.42_count<11, "Partners 0-10", "Partners 11 or more")

# Make a heatmap with key evidence
# Select key variables and reformat for plotting
polyq_genes_info_heatmap <- polyq_genes_info %>%
  dplyr::select(gene, interacting_partners, 
                rna_tissue_specificity, pLI, druggable_any)

# Melt
polyq_genes_info_heatmap  <- reshape2::melt(polyq_genes_info_heatmap, id="gene") %>%
  dplyr::rename(feature = variable)

# Perform replacements to group under favorable, unfavorable, etc.
polyq_genes_info_heatmap$value <- gsub("Tissue enhanced", 1, polyq_genes_info_heatmap$value)
polyq_genes_info_heatmap$value <- gsub("Tissue enriched", 1, polyq_genes_info_heatmap$value)
polyq_genes_info_heatmap$value <- gsub("Low tissue specificity", 0, polyq_genes_info_heatmap$value)
polyq_genes_info_heatmap$value <- gsub("Intolerant", 0, polyq_genes_info_heatmap$value)
polyq_genes_info_heatmap$value <- gsub("Tolerant", 1, polyq_genes_info_heatmap$value)
polyq_genes_info_heatmap$value <- gsub("No", 0, polyq_genes_info_heatmap$value)
polyq_genes_info_heatmap$value <- gsub("Yes", 1, polyq_genes_info_heatmap$value)
polyq_genes_info_heatmap$value <- gsub("Partners 0-10", 1, polyq_genes_info_heatmap$value)
polyq_genes_info_heatmap$value <- gsub("Partners 11 or more", 0, polyq_genes_info_heatmap$value)
polyq_genes_info_heatmap$value <- gsub("Partners 0-10", 1, polyq_genes_info_heatmap$value)
polyq_genes_info_heatmap$value <- gsub("Partners 11 or more", 0, polyq_genes_info_heatmap$value)

# Replace AR pli with 0.5 score for later plot ordering since it unknown
polyq_genes_info_heatmap[21, 3] <- 0.5

# Add a plot value
polyq_genes_info_heatmap$plot_value <- gsub(1, "Favorable", polyq_genes_info_heatmap$value)
polyq_genes_info_heatmap$plot_value <- gsub(0.5, "Unresolved", polyq_genes_info_heatmap$plot_value)
polyq_genes_info_heatmap$plot_value <- gsub(0, "Unfavorable", polyq_genes_info_heatmap$plot_value)

# Re-level feature factor for plotting
gene_totals <- polyq_genes_info_heatmap %>%
  group_by(gene) %>%
  summarise(total_value = sum(as.numeric(as.character(value)), na.rm = TRUE))

desired_order <- c( "interacting_partners", "rna_tissue_specificity",
                   "pLI", "druggable_any")


# Re-level the feature variable with the desired order
polyq_genes_info_heatmap$feature <- factor(polyq_genes_info_heatmap$feature, levels = desired_order)
gene_totals <- gene_totals[order(gene_totals$total_value, decreasing = TRUE), ]
polyq_genes_info_heatmap$gene <- factor(polyq_genes_info_heatmap$gene, levels = gene_totals$gene)
rm(desired_order)
rm(gene_totals)

# Make a heatmap
polyq_genes_info_heatmap %>%
  ggplot(aes(x=reorder(gene, value==0), y=(feature), fill=plot_value)) + 
  geom_tile(color="white", linewidth=0.5) +
  coord_equal() +
  theme_minimal() +
  labs(x=NULL, y=NULL) +
  scale_fill_manual(
    values = c("#21918c", "#440154", "grey"),
    labels = c("Favorable", "Unfavorable", "Unresolved")
  ) +theme(legend.position="top", legend.title=element_blank(), 
           axis.text.y = element_text(size=12, face="bold"),
           axis.text.x = element_text(size=12, face="bold.italic")) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  scale_y_discrete(labels = c("Interacting Partners", "RNA tissue specificity", "pLi tolerance", "Druggable any"))

ggsave("/research/2023_polyQ/results/polyq_genes_druggability.pdf", height=5.36, width=7.58*0.8, units='in')


