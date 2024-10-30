# Project
This contains the R and BASH scripts, results and data required for the analyses performed for **"Unbiased genomic characterization of human traits associated with polyglutamine disease genes to inform therapeutic strategies."**

## Directory Structure
-`data/`: Contains data used for analyses. Some data sets are restricted and have not been included in the folder. 
  -`categories_dgidb_2024_07_23.tsv`: Data for the drug categories of various genes exported from the [Drug Gene Interaction Database (DGIdb)](https://www.dgidb.org/) on July 23, 2024.
  - `gnomad.v2.1.1.lof_metrics.by_gene.txt`: *Not included.* Loss-of-function metrics by gene from gnomAD v2.1.1. Available from: [Download](https://gnomad.broadinstitute.org/downloads)
  -`polyQ_otg_29_09-23.tsv`: Data extracted from the Open Targets Genetic portal, specifically related to the polyglutamine disease genes (accessed January 24, 2024).
  -`polyQ_genes_molecular_interactions_interactors.txt`: IntAct molecular interaction data extracted from the Open Targets database for polyglutamine disease genes (Data extracted on October 25, 2024).
  - `proteinatlas.tsv`: *Not included.* Data from the Human Protein Atlas. Available from: [Download](https://www.proteinatlas.org/about/download)

-`results/`: Contains figures generated from the analyses for the current study. 

-`scr/`: Contains scripts used in the analysis and to generate figures.
 -`R/`
  -`1_polyq_analyses.R`: code to curate the polyglutamine gene-associated traits from the Open Targets Genetics database, perform gene-related analyses and generate figures. 
  -`2_polyq_druggability_analyses.R`: code to perform polyglutamine gene druggability-related analysis and generate figures.
  -`3_polyq_annotation_analyses.R`: code to generate data set used to retrieve the annotations of the variants associated with the polyglutamine trait associations.
 -`bash/`
  -`1_polyq_otg_annotation.sh`: code to retrieve annotations of the variants associated with the polyglutamine gene-trait associations using Ensemble Variant Effect predictor.
