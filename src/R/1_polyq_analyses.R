# Install relevant libraries
install.packages("ghql")
install.packages("ggrepel")
install.packages("BiocManager")
BiocManager::install("biomaRt") 
install.packages("dplyr")
install.packages("cli")

# Load the libraries 
library(dplyr)
library(data.table) 
library(ghql)
library(jsonlite) 
library(ggplot2)
library(ggrepel)
library(readxl)
library(biomaRt) 
library(forcats)
library(viridisLite)
library(viridis)
library(patchwork)
library(reshape)
library(arcdiagram)
library(igraph)
library(ggraph)

# These scripts were adapted from:
# https://blog.opentargets.org/summer-school-crash-course-part-2/
# https://community.opentargets.org/t/accessing-locus-to-gene-l2g-and-colocalisation-data-from-the-genetics-portal-by-querying-the-api-using-python/124
# See schema for what you can extract
# https://api.genetics.opentargets.org/graphql/schema
# Set up to query Open Targets Genetics API
otg_cli <- GraphqlClient$new(url = "https://api.genetics.opentargets.org/graphql")
otg_qry <- Query$new()

# Query for GWAS study locus details
otg_qry$query('l2g_query', 'query GenePageQuery($geneId: String!){
    geneInfo(geneId: $geneId){
    id
    symbol
  }
   studiesAndLeadVariantsForGeneByL2G(geneId: $geneId) {
    pval
    yProbaModel
    study {
      studyId
      traitReported
      traitCategory
      pubAuthor
      pubDate
      pmid
      nInitial
      nReplication
      hasSumstats
      numAssocLoci
    }
    variant {
      chromosome
      position
      refAllele
      altAllele
      rsId
      id
    }
    odds {
      oddsCI
      oddsCILower
      oddsCIUpper
    }
    beta {
      betaCI
      betaCILower
      betaCIUpper
      direction
    }
  }
}')

## Execute the query 
# Write a function that fetches L2Gs for the polyQ genes
# Use this later to automate as you can only fetch one gene at a time
fetch_l2g <- function(current_geneId) {
  variables = list(geneId = current_geneId)
  result_gene <- fromJSON(otg_cli$exec(otg_qry$queries$l2g_query, variables, flatten = TRUE))$data
  if(nrow(as.data.frame(result_gene$studiesAndLeadVariantsForGeneByL2G))==0)
    result_gene_df <- NULL
  else
    result_gene_df <- as.data.frame(result_gene) %>% flatten()
  return(result_gene_df)
}
colnames(result_df)

# Load the polyQ gene information
polyQ_disorder_genes <- read_excel("/research/2023_polyQ/otg/data/polyQ_disorders_2023.xlsx",sheet = "Sheet1")

# Retreive the ensembl ids for the polyQ genes
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", archive=FALSE, verbose=TRUE)
filters <- listFilters(ensembl)
attributes <- listAttributes(ensembl)
filterlist <- unique(polyQ_disorder_genes$Gene) # unique function returns a vector or dataframe without duplicated rows
polyq_gene_info <- getBM(attributes = c("hgnc_symbol","ensembl_gene_id", "start_position", "end_position"),filters = c("hgnc_symbol"),values = filterlist, mart = ensembl)

# Extract L2Gs for the polyQ genes 
query_genes <- polyq_gene_info$ensembl_gene_id
l2g_combined <- NULL
for (gene in query_genes) {
  print(gene)
  temp_l2g <- fetch_l2g(gene)
  l2g_combined <- rbind(l2g_combined, temp_l2g) 
}

# Clean and rename the column names
colnames(polyq_gene_info)
colnames(l2g_combined)
colnames(l2g_combined) <- gsub("studiesAndLeadVariantsForGeneByL2G.", "", colnames(l2g_combined))
colnames(l2g_combined) <- gsub("geneInfo.", "", colnames(l2g_combined))
colnames(l2g_combined) <- gsub("yProbaModel", "L2G", colnames(l2g_combined))

# Save the polyQ otg data
# Data was downloaded on 29th September 2023
write.table(l2g_combined, file = "/research/2023_polyQ/otg/data/polyQ_otg_29_09-23.tsv", sep = "\t", row.names = FALSE)
l2g_combined <- fread("/research/2023_polyQ/otg/data/polyQ_otg_29_09-23.tsv")

# Filter for select associations (remove count based, etc)
l2g_combined_select <- l2g_combined %>% 
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

# Filter for the L2G significant (L2G>0.5) associations 
l2g_combined_sig <- l2g_combined %>%
  filter(L2G>0.5)
l2g_combined_select_sig <- l2g_combined_select %>%
  filter(L2G>0.5)

## Plot the polyQ L2G significant associations
l2g_combined_sig %>%
  ggplot(aes(x = fct_infreq(symbol,), fill = symbol)) +
  theme_minimal() + theme_bw() +
  geom_bar(width = 0.6) +
  scale_fill_manual(values=c("THAP11" = "#addc30","TBP" ="#fde725","ATN1" = "#5ec962", "ATXN7" = "#28ae80", "ATXN3" = "#21918c","ATXN1" = "#2c728e","CACNA1A" = "#3b528b","HTT" = "#472d7b","ATXN2" = "#440154"))+
  labs(x = "Gene",
       y ="Number of traits reported") +
  geom_text(aes(label = ..count..), stat = "count", vjust = -0.2, colour = "black", size = 3.5) +
  #labs(title = "Total PolyQ L2G Studies")+
  theme(plot.title = element_text(hjust = 0.5, size = 9),
        axis.text.x = element_text(face = "bold.italic", angle = 90, size = 12),
        axis.text.y = element_text(face = "bold", size = 12),
        legend.position = "none",
        axis.title.x = element_text(size = 10, face="bold"),
        axis.title.y = element_text(size = 10, face="bold")) +
  coord_cartesian(ylim = c(0, 75))

ggsave("/research/2023_polyQ/otg/results/polyq_genes_l2g_significant_associations.pdf", height=5.36, width=5.36, units='in')

## Plot the measurement category vs other categories 
measure <- l2g_combined_sig %>%
  mutate(measure = ifelse(study.traitCategory == "measurement", "measurement", "non-measurement"))
measure %>%
  ggplot(aes(x = fct_infreq(symbol,), fill = measure)) +
  theme_minimal() + theme_bw() +
  geom_bar(width = 0.6) +
  scale_fill_manual(values=c("measurement" = "#63d6d2", "non-measurement" = "#21111b"))+
  labs(x = "Gene",
       y ="Number of traits reported") +
  theme(axis.text.x = element_text(angle = 90,  size = 12 ,face = "bold.italic"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 10, face = "bold"),
        legend.position = "bottom",
        legend.title = element_blank())

ggsave("/research/2023_polyQ/otg/results/polyq_genes_l2g_sig_measure_propotion.pdf", height=5.36, width=5.36, units='in')

## Plot an L2G heatmap showing distinct non-measurement categories for each gene
heatmap <- l2g_combined_sig %>%
  dplyr::select(c(symbol, study.traitCategory)) %>%
  filter(!(study.traitCategory %like% "measurement")) %>%
  group_by(symbol) %>%
  arrange(study.traitCategory) %>%
  add_count(study.traitCategory, name = "total_trait_categories")

# Calculate the total number of genes for each category
category_gene_counts <- heatmap %>%
  group_by(study.traitCategory) %>%
  summarise(total_genes = sum(total_trait_categories))

# Reorder the levels of the study.traitCategory factor variable based on the total number of genes for each category
heatmap$study.traitCategory <- factor(heatmap$study.traitCategory, 
                                      levels = category_gene_counts$study.traitCategory[order(-category_gene_counts$total_genes)])

# Define the order of genes
gene_order <- c("ATXN1", "HTT", "ATXN2", "ATXN7", "CACNA1A")

# Reorder the levels of the symbol factor variable according to the specified gene order
heatmap$symbol <- factor(heatmap$symbol, levels = gene_order)

## Plot the heatmap with rearranged genes and categories
ggplot(heatmap, aes(x = symbol, y = study.traitCategory, fill = total_trait_categories)) +
  geom_tile(color = "white",
            lwd = 0.5,
            linetype = 1) +
  theme_bw() +
  geom_text(aes(label = total_trait_categories), color = "white", size = 4) +
  theme(axis.text.x = element_text(angle = 90, size = 12, face = "bold.italic"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 10, face = "bold"),
        legend.position = "none") +
  xlab("Gene") + 
  ylab("Study traits per Category") +
  labs(fill = "Total Categories") +
  coord_fixed()

ggsave("/research/2023_polyQ/otg/results/polyq_genes_l2g_sig_heatmap.pdf", height=4.36, width=5.36, units='in')

## Plot the network graph
# Create the dataframe to use for the network
l2g_combined_network <- l2g_combined_sig %>%
  dplyr::select(c("symbol", "study.traitReported", "study.traitCategory"))
l2g_combined_network1 <- l2g_combined_sig %>%
  dplyr::select(c("symbol", "study.traitReported", "study.traitCategory"))

# Create the edge list 
edge_df <- inner_join(l2g_combined_network, l2g_combined_network1, by = "study.traitReported") 

# Remove all the duplicate rows
edge_df <- edge_df %>% dplyr::filter(symbol.x != symbol.y)

# Calculate frequency of genes
gene_counts <- sort(table(edge_df$symbol.x), decreasing = TRUE)

# Sort genes to the desired order
sorted_genes <- names(gene_counts)

# Reorder edge_df and select the desired columns
edge_df <- edge_df[order(factor(edge_df$symbol.x, levels = sorted_genes), 
                         factor(edge_df$symbol.y, levels = sorted_genes)), ]
edge_df <- edge_df %>%
  dplyr::select(c(symbol.x, symbol.y))

# Convert to matrix
edge_df <- as.matrix(edge_df)

# Define specific colors for the nodes
node_colors <- c("ATXN1" = "#2c728e", "HTT" = "#472d7b", "ATXN2" = "#440154", 
                 "ATXN7" = "#28ae80", "ATXN3" = "#21918c")

# Calculate node degree (number of associations) for node sizing 
node_degrees <- table(c(edge_df))
node_degrees <- sort(node_degrees, decreasing = TRUE)

# Scale degrees to a range between 1 and 
scaled_degrees <- (node_degrees - min(node_degrees)) / (max(node_degrees) - min(node_degrees)) * 4 + 1

# Map the defined colors
sorted_nodes <- names(sort(node_degrees, decreasing = TRUE))
valid_nodes <- sorted_nodes[sorted_nodes %in% names(node_colors)]
mapped_colors <- node_colors[valid_nodes]

# Calculate line widths for arcs based on the number of associations
# (More associations = thicker arc)
edge_weights <- apply(edge_df, 1, function(x) sum(x[1] %in% sorted_nodes) + sum(x[2] %in% sorted_nodes))
arc_widths <- (edge_weights - min(edge_weights)) / (max(edge_weights) - min(edge_weights)) * 2 + 1

# Open a PDF device
pdf("/research/2023_polyQ/otg/results/polyq_genes_l2g_sig_network.pdf", height=4.86, width=6.00) 

# Plot the network graph
arcplot(edge_df, col.arcs = hsv(0, 0, 0.2, 0.25),
        lwd.arcs = arc_widths, col.nodes = mapped_colors,
        cex.nodes = sqrt(scaled_degrees[valid_nodes]),
        cex.labels = 0.7, font = 3, srt = 45, line = 0.8)

# Close the PDF device
dev.off()

## Plot the polyQ select and significat associations 
ggplot(l2g_combined_select_sig, aes(fct_infreq(symbol,), fill = symbol)) +
  geom_bar(width = 0.4) + theme_minimal() + theme_bw() +
  scale_fill_manual(values=c("ATXN7" = "#28ae80","ATXN1" = "#2c728e","CACNA1A" = "#3b528b","HTT" = "#472d7b","ATXN2" = "#440154"))+
  labs(x = "Gene", y = "Number of traits reported") +
  geom_text(aes(label = ..count..), stat = "count", vjust = -0.2, colour = "black", size = 3.5) +
  theme(plot.title = element_text(hjust = 0.5, size = 9),
        axis.text.x = element_text(face = "bold.italic", angle = 90, size = 12),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10),
        legend.position = "none") 

ggsave("/research/2023_polyQ/otg/results/polyq_genes_l2g_high_confidence_associations.pdf", height=5.36, width=5.36, units='in')

## Plot polyQ Lollipop for the significant associations 
# Reformat and take out the excess wording from the traits 
current <- l2g_combined_select_sig
current$study.traitReported <- gsub(" \\[EA\\])", "", current$study.traitReported)
current$study.traitReported <- gsub(" \\(years of education\\)", "", current$study.traitReported)
current$study.traitReported <- gsub(" \\(special factor of neuroticism\\)", "", current$study.traitReported)
current$study.traitReported <- gsub(" \\(healthspan, parental lifespan or longevity)", "", current$study.traitReported)
current$study.traitReported <- gsub(" \\(multivariate analysis\\)", "", current$study.traitReported)
current$study.traitReported <- gsub(" \\[MTAG\\]", "", current$study.traitReported)
current$study.traitReported <- gsub(" \\(MTAG\\)", "", current$study.traitReported)
current$study.traitReported <- gsub(" \\(Hounsfield unit scale\\)", "", current$study.traitReported)
current$study.traitReported <- gsub(" \\(pleiotropy\\)", "", current$study.traitReported)
current$study.traitReported <- gsub(" \\(years of education\\)", "", current$study.traitReported)
current$study.traitReported <- gsub(" \\[EA\\]", "", current$study.traitReported)
current$study.traitReported <- gsub(" \\[conditional-joint\\]", "", current$study.traitReported)
current$study.traitReported <- gsub(" \\[MTAG\\]", "", current$study.traitReported)
current$study.traitReported <- gsub(" \\(MTAG\\)", "", current$study.traitReported)
current$study.traitReported <- gsub("/atrial flutter", "", current$study.traitReported)

# Filter out the duplicate traits and restructure the data
current <- current %>%
  group_by(study.traitReported) %>%
  slice_max(order_by = L2G) %>%
  ungroup()

# Some duplicate traits have the same L2G so take out the one that appears first
current <- current %>%
  group_by(study.traitReported) %>%
  slice_head(n = 1) %>%
  ungroup()

# Plot the lollipop plot
current %>%
  group_by(symbol) %>%
  arrange(L2G) %>%   
  mutate(study.traitReported = factor (study.traitReported, levels = study.traitReported)) %>%   # This trick update the factor levels
  ggplot( aes(study.traitReported, L2G, color = symbol)) +
  theme_minimal() + theme_bw() +
  geom_segment( aes(x = study.traitReported, xend = study.traitReported, y = 0, yend = L2G)) +
  geom_point( size = 2.6) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", size = 0.4) +
  scale_color_manual(values = c("ATXN7" = "#28ae80","ATXN1" = "#2c728e","CACNA1A" = "#3b528b","HTT" = "#472d7b","ATXN2" = "#440154")) +
  coord_flip() +
  labs(x = "Study traits reported") +
  labs(color = "Gene") + 
  theme(axis.title = element_text(face = "bold"),
        title = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "right")

ggsave("/research/2023_polyQ/otg/results/polyq_genes_l2g_significant_lollipop.pdf", height=5.36, width=8.36, units='in')
