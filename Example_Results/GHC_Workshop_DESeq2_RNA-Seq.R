# DESeq2 Analysis for Plant Regen RNA-Seq Experiment======


# Dependencies ----
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(knitr)
library(gprofiler2)
library(tidyverse)
library(gridExtra)
library(htmltools)
library(DESeq2)
library(ashr)

# Visualization Dependencies ----
library(RColorBrewer)
library(gplots)
library(SummarizedExperiment)
library(rtracklayer)
library(plotly)

# GO analysis Dependencies ----
library(UpSetR)
library(clusterProfiler)
library(pathview)
library(enrichplot)
library(ggplot2)

# Environment Cleanup----------

rm(list = ls())  # Ensures no contamination from previous sessions


# Working Directory----------

setwd("~/Desktop/Research/Dissertation_Work/Chapter3_Regeneration/Regen_Trial_Dir")


# Load Annotation Data---------

# TAIR annotation provides gene names and functional metadata
ATL_Annotation <- read.csv(
  "ATL_v3.functional_annotation.csv",
  sep = ',',
  header = TRUE
)


# Load RNA-Seq Count Matrix------

# Rows = genes, Columns = samples
RNA_Counts_df <- as.data.frame(
  read.csv("RegenTrial_Gene_Count_Matrix.txt", sep = "\t", header = TRUE)
)

head(RNA_Counts_df)

# Assign biologically meaningful column names
colnames(RNA_Counts_df) 

# Prepare Count Matrix--------

RNA_Counts_mat <- RNA_Counts_df
rownames(RNA_Counts_mat) <- RNA_Counts_mat$Geneid
RNA_Counts_mat$Geneid <- NULL

# Remove genes with zero total counts
geneTotals <- rowSums(RNA_Counts_mat)
countsNonZero <- RNA_Counts_mat[geneTotals > 0, ]

# Remove genes with missing values
countsNonZero <- na.omit(countsNonZero)

# Remaining genes after filtering
nrow(countsNonZero)


# Load Metadata---------

metaData <- read.csv("GHC_Regen_Metadata.csv")

colnames(metaData)

# Inspect factor levels (critical for DESeq2 reference level interpretation)
factor(metaData$Plate)
factor(metaData$Timepoint)


# Construct DESeq2 Dataset-------

dds <- DESeqDataSetFromMatrix(
  countsNonZero,
  colData = metaData,
  design = ~Timepoint
)

dds$Timepoint <- factor(dds$Timepoint, levels = c("1D","8D","14D","28D"))
dds$Timepoint <- relevel(dds$Timepoint, ref ="1D")
model.matrix(design(dds), colData(dds))

# Design matrix inspection
dds_model <- as.data.frame(model.matrix(design(dds), colData(dds)))

rownames(dds_model) <- colnames(RNA_Counts_mat)


# Differential Expression Analysis------

dds <- DESeq(dds)

# Inspect coefficient names available for contrasts
resultsNames(dds)


# Extract Contrast: 8D vs 1D------

D1_vs_D8 <- results(dds, name = "Timepoint_8D_vs_1D")
summary(D1_vs_D8)

# Explicit contrast definition (preferred for clarity)
D1_vs_D8 <- results(dds, contrast = c("Timepoint", "8D", "1D"))
summary(D1_vs_D8)

#Metadata describing result columns
mcols(D1_vs_D8, use.names = TRUE)


#### Convert Results to DataFrame-----
D1_vs_D8_df <- as.data.frame(
  results(dds, contrast = c("Timepoint", "8D", "1D"))
)

D1_vs_D8_df$Gene_ID <- row.names(D1_vs_D8_df)

# Ensure Gene_ID is first column
D1_vs_D8_df <- D1_vs_D8_df[, c(
  "Gene_ID",
  setdiff(colnames(D1_vs_D8_df), "Gene_ID")
)]


#### Filter Significant Genes--------

# Criteria: adjusted p-value < 0.05 and |log2FC| > 2
Sig_genes_D1_vs_D8_df <- as.data.frame(
  D1_vs_D8_df[
    D1_vs_D8_df$padj < 0.05 &
      abs(D1_vs_D8_df$log2FoldChange) > 2,
  ]
)

Sig_genes_D1_vs_D8_df <- na.omit(Sig_genes_D1_vs_D8_df)


#### Merge DEGs with Annotation file--------

# Attach biological meaning (gene names, functions)
D1_vs_D8_Annotated_df <- merge(
  ATL_Annotation,
  D1_vs_D8_df,
  by.x = "Gene",
  by.y = "Gene_ID"
)

D1_vs_D8_SIG_Annotated_df <- merge(
  ATL_Annotation,
  Sig_genes_D1_vs_D8_df,
  by.x = "Gene",
  by.y = "Gene_ID"
)

# Remove redundant merge column
D1_vs_D8_Annotated_df <- D1_vs_D8_Annotated_df[
  , !names(D1_vs_D8_Annotated_df) %in% c("Gene_ID_raw")
]

D1_vs_D8_SIG_Annotated_df <- D1_vs_D8_SIG_Annotated_df[
  , !names(D1_vs_D8_SIG_Annotated_df) %in% c("Gene_ID_raw")
]


#### Sort DEGs by Effect Size (Log2FC)---------

D1_vs_D8_Annotated_df <- D1_vs_D8_Annotated_df[
  order(-D1_vs_D8_Annotated_df$log2FoldChange),
]

D1_vs_D8_SIG_Annotated_df <- D1_vs_D8_SIG_Annotated_df[
  order(-D1_vs_D8_SIG_Annotated_df$log2FoldChange),
]


#### Summary Counts 8D vs 1D------

Unique_All_Genes_Amount_D1_vs_D8 <- length(
  unique(D1_vs_D8_Annotated_df$Gene)
)

Unique_Sig_Genes_Amount_D1_vs_D8 <- length(
  unique(D1_vs_D8_SIG_Annotated_df$Gene)
)


Unique_All_Genes_Amount_D1_vs_D8
Unique_Sig_Genes_Amount_D1_vs_D8



# Extract Contrast: Timepoint 14D vs 1D ------

D1_vs_D14 <- results(dds, name = "Timepoint_14D_vs_1D")
summary(D1_vs_D14)

# Explicit contrast definition (preferred for clarity)
D1_vs_D14 <- results(dds, contrast = c("Timepoint", "14D", "1D"))
summary(D1_vs_D14)

# Metadata describing result columns
mcols(D1_vs_D14, use.names = TRUE)


#### Convert Results to DataFrame-------

D1_vs_D14_df <- as.data.frame(
  results(dds, contrast = c("Timepoint", "14D", "1D"))
)

D1_vs_D14_df$Gene_ID <- row.names(D1_vs_D14_df)

# Ensure Gene_ID is first column
D1_vs_D14_df <- D1_vs_D14_df[, c(
  "Gene_ID",
  setdiff(colnames(D1_vs_D14_df), "Gene_ID")
)]


#### Filter Significant Genes--------

# Criteria: adjusted p-value < 0.05 and |log2FC| > 0.5
Sig_genes_D1_vs_D14_df <- as.data.frame(
  D1_vs_D14_df[
    D1_vs_D14_df$padj < 0.05 &
      abs(D1_vs_D14_df$log2FoldChange) > 2,
  ]
)

Sig_genes_D1_vs_D14_df <- na.omit(Sig_genes_D1_vs_D14_df)


#### Merge DEGs with Annotation file--------

# Attach biological meaning (gene names, functions)
D1_vs_D14_Annotated_df <- merge(
  ATL_Annotation,
  D1_vs_D14_df,
  by.x = "Gene",
  by.y = "Gene_ID"
)

D1_vs_D14_SIG_Annotated_df <- merge(
  ATL_Annotation,
  Sig_genes_D1_vs_D14_df,
  by.x = "Gene",
  by.y = "Gene_ID"
)

# Remove redundant merge column
D1_vs_D14_Annotated_df <- D1_vs_D14_Annotated_df[
  , !names(D1_vs_D14_Annotated_df) %in% c("Gene_ID_raw")
]

D1_vs_D14_SIG_Annotated_df <- D1_vs_D14_SIG_Annotated_df[
  , !names(D1_vs_D14_SIG_Annotated_df) %in% c("Gene_ID_raw")
]


#### Sort DEGs by Effect Size (Log2FC)---------

D1_vs_D14_Annotated_df <- D1_vs_D14_Annotated_df[
  order(-D1_vs_D14_Annotated_df$log2FoldChange),
]

D1_vs_D14_SIG_Annotated_df <- D1_vs_D14_SIG_Annotated_df[
  order(-D1_vs_D14_SIG_Annotated_df$log2FoldChange),
]


#### Summary Count 14D vs 1D------

Unique_All_Genes_Amount_D1_vs_D14 <- length(
  unique(D1_vs_D14_Annotated_df$Gene)
)

Unique_Sig_Genes_Amount_D1_vs_D14 <- length(
  unique(D1_vs_D14_SIG_Annotated_df$Gene)
)


Unique_All_Genes_Amount_D1_vs_D14
Unique_Sig_Genes_Amount_D1_vs_D14



# Extract Contrast: Timepoint 28D vs 1D ------

D1_vs_D28 <- results(dds, name = "Timepoint_28D_vs_1D")
summary(D1_vs_D28)

# Explicit contrast definition (preferred for clarity)
D1_vs_D28 <- results(dds, contrast = c("Timepoint", "28D", "1D"))
summary(D1_vs_D28)

# Metadata describing result columns
mcols(D1_vs_D28, use.names = TRUE)


#### Convert Results to DataFrame-------

D1_vs_D28_df <- as.data.frame(
  results(dds, contrast = c("Timepoint", "28D", "1D"))
)

D1_vs_D28_df$Gene_ID <- row.names(D1_vs_D28_df)

# Ensure Gene_ID is first column
D1_vs_D28_df <- D1_vs_D28_df[, c(
  "Gene_ID",
  setdiff(colnames(D1_vs_D28_df), "Gene_ID")
)]


#### Filter Significant Genes--------

# Criteria: adjusted p-value < 0.05 and |log2FC| > 0.5
Sig_genes_D1_vs_D28_df <- as.data.frame(
  D1_vs_D28_df[
    D1_vs_D28_df$padj < 0.05 &
      abs(D1_vs_D28_df$log2FoldChange) > 2,
  ]
)

Sig_genes_D1_vs_D28_df <- na.omit(Sig_genes_D1_vs_D28_df)


#### Merge DEGs with Annotation file--------

# Attach biological meaning (gene names, functions)
D1_vs_D28_Annotated_df <- merge(
  ATL_Annotation,
  D1_vs_D28_df,
  by.x = "Gene",
  by.y = "Gene_ID"
)

D1_vs_D28_SIG_Annotated_df <- merge(
  ATL_Annotation,
  Sig_genes_D1_vs_D28_df,
  by.x = "Gene",
  by.y = "Gene_ID"
)

# Remove redundant merge column
D1_vs_D28_Annotated_df <- D1_vs_D28_Annotated_df[
  , !names(D1_vs_D28_Annotated_df) %in% c("Gene_ID_raw")
]

D1_vs_D28_SIG_Annotated_df <- D1_vs_D28_SIG_Annotated_df[
  , !names(D1_vs_D28_SIG_Annotated_df) %in% c("Gene_ID_raw")
]


#### Sort DEGs by Effect Size (Log2FC)---------

D1_vs_D28_Annotated_df <- D1_vs_D28_Annotated_df[
  order(-D1_vs_D28_Annotated_df$log2FoldChange),
]

D1_vs_D28_SIG_Annotated_df <- D1_vs_D28_SIG_Annotated_df[
  order(-D1_vs_D28_SIG_Annotated_df$log2FoldChange),
]


#### Summary Counts 28D vs 1D------

Unique_All_Genes_Amount_D1_vs_D28 <- length(
  unique(D1_vs_D28_Annotated_df$Gene)
)

Unique_Sig_Genes_Amount_D1_vs_D28 <- length(
  unique(D1_vs_D28_SIG_Annotated_df$Gene)
)


Unique_All_Genes_Amount_D1_vs_D28
Unique_Sig_Genes_Amount_D1_vs_D28


# Filter to representative gene IDs (remove redundant isoforms) ----

Atl_Repr_df <- read.csv("GHC_Representative_Atlv3-Ids.txt", header = TRUE, sep = '\t'
)

# Keep only genes present in representative set (ensures 1 gene model per locus)
D1_vs_D8_SIG_rep  <- inner_join(D1_vs_D8_SIG_Annotated_df,  Atl_Repr_df, by = "Gene")
D1_vs_D14_SIG_rep <- inner_join(D1_vs_D14_SIG_Annotated_df, Atl_Repr_df, by = "Gene")
D1_vs_D28_SIG_rep <- inner_join(D1_vs_D28_SIG_Annotated_df, Atl_Repr_df, by = "Gene")


# Count unique representative genes per comparison ----

Unique_Repr_Sig_Genes_Amount_D1_vs_D8  <- length(unique(D1_vs_D8_SIG_rep$Gene))
Unique_Repr_Sig_Genes_Amount_D1_vs_D14 <- length(unique(D1_vs_D14_SIG_rep$Gene))
Unique_Repr_Sig_Genes_Amount_D1_vs_D28 <- length(unique(D1_vs_D28_SIG_rep$Gene))

# Compare original vs representative counts (sanity check for filtering impact)
Unique_Sig_Genes_Amount_D1_vs_D8
Unique_Repr_Sig_Genes_Amount_D1_vs_D8

Unique_Sig_Genes_Amount_D1_vs_D14
Unique_Repr_Sig_Genes_Amount_D1_vs_D14

Unique_Sig_Genes_Amount_D1_vs_D28
Unique_Repr_Sig_Genes_Amount_D1_vs_D28


# Combine all DEGs into single table with comparison label ----

Regen_DEGs <- bind_rows(
  D1_vs_D8_SIG_rep  %>% mutate(Comparison = "D1_vs_D8"),
  D1_vs_D14_SIG_rep %>% mutate(Comparison = "D1_vs_D14"),
  D1_vs_D28_SIG_rep %>% mutate(Comparison = "D1_vs_D28")
)




# Basic Visualization: PCA Analysis using DESeq2 VST--------


# Metadata sanity check-------

# Ensure metadata is complete and corresponds to all samples used in analysis
metaData 



# Variance Stabilizing Transformation (VST)------

# Converts count data into a homoscedastic space (variance stabilized)
# Suitable for PCA, clustering, and visualization (NOT for DE testing)
dds_vst <- vst(dds, blind = FALSE)

# Inspect available sample annotations for grouping variables
colData(dds_vst)



# Align Metadata with Expression Data------

# Ensures correct mapping between sample metadata and transformed counts

# Verify required metadata fields exist
stopifnot("Library" %in% names(metaData))


# Reorder metadata to match column order of expression matrix
PCA_metaData <- metaData[
  match(colnames(dds_vst), metaData$Library),
]

# Assign rownames for consistent indexing
rownames(PCA_metaData) <- PCA_metaData$Library

# Attach Metadata to DESeq Object-------
# Required so plotPCA can use grouping variables
colData(dds_vst)$Plate  <- PCA_metaData$Plate
colData(dds_vst)$Timepoint <- PCA_metaData$Timepoint



# PCA Computation---------

# Uses top variable genes (ntop) to compute principal components
pcaData <- plotPCA(
  dds_vst,
  intgroup = c("Timepoint"),
  returnData = TRUE,
  ntop = 2500
)

# Factor Level Ordering---------
# Controls plotting order and legend appearance

# Ensure correct ordering
pcaData$Timepoint <- factor(
  pcaData$Timepoint,
  levels = c("1D","8D","14D","28D")
)

pcaData$Plate <- factor(pcaData$Plate)

# Extract percentage of variance explained by PCs
percentVar <- round(100 * attr(pcaData, "percentVar"))



# Visualization Setup------

# Custom color palette for genotype groups
color_palette <- c("#135E4B", "#4CB572","#A1D8B5","#CCDCDB")



# PCA Plot--------

# PCA plot: visualize sample clustering from global expression (PC1/PC2 = main variance axes)

ggplot(pcaData, aes(PC1, PC2, color = Timepoint, shape = Plate)) +
  
  geom_point(size = 3, stroke = 0.8) +
  
  scale_color_manual(
    name = "Timepoint",
    values = color_palette,
    breaks = c("1D","8D","14D","28D")  # legend order
  ) +
  
  scale_shape_manual(
    name = "Plate",
    values = c("CIM" = 16, "SIM" = 17)  # circle vs triangle
  ) +
  
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  
  coord_fixed() +
  
  theme_light() +
  
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title = element_text(size = 12),
    legend.position = "right"
  ) +
  
  ggtitle("PCA Analysis")







# Load TAIR annotation (for gene names + functional info) ----

TAIR_Annot <- read.csv("GHC_TAIR10_Annotation_Dataframe.tsv", sep = '\t')


# Load GO annotation (Soltu → GO + TAIR mapping) ----

ATL_GO <- read.csv(
  "ATL_v3.working_models.go_slim.obo",
  sep = '\t',
  header = FALSE
)


# Clean and standardize GO annotation ----

colnames(ATL_GO) <- c(
  "Source",
  "Gene",              # Soltu gene ID
  "Gene_version",
  "NA_col",
  "GO_ID",
  "TAIR_raw",
  "Evidence",
  "NA1","Aspect","NA2","NA3","Type","Taxon","Date","Source2"
)

ATL_GO <- ATL_GO %>%
  mutate(
    TAIR_ID = str_remove(TAIR_raw, "^TAIR:"),   # normalize TAIR IDs for joining
    GO_ID   = str_replace(GO_ID, "\\.", ":")    # fix GO format consistency
  ) %>%
  select(Gene, GO_ID, TAIR_ID)


# Merge GO + TAIR annotation into DEG table ----

Regen_DEGs_TAIR_Annot <- Regen_DEGs %>%
  left_join(ATL_GO, by = "Gene", relationship = "many-to-many")  # 1 gene → many GO terms

Regen_DEGs_TAIR_Annot <- Regen_DEGs_TAIR_Annot %>%
  left_join(
    TAIR_Annot %>%
      select(Gene_ID, Gene_name, TAIR_Function = Functional_Annotation),
    by = c("TAIR_ID" = "Gene_ID")
  )

#Annotating the whole RNA_Counts_df for future references
RNA_counts_rep<- RNA_Counts_df
RNA_counts_annot <- RNA_counts_rep %>%
  left_join(ATL_GO, by =c("Geneid" = "Gene"), relationship = "many-to-many") %>%
  left_join(
    TAIR_Annot %>%
      select(Gene_ID, Gene_name, TAIR_Function = Functional_Annotation),
    by = c("TAIR_ID" = "Gene_ID")
  )

ATL_GO_collapsed <- ATL_GO %>%
  group_by(Gene) %>%
  summarise(
    GO_ID = paste(unique(GO_ID), collapse = "/"),
    TAIR_ID = first(TAIR_ID),
    .groups = "drop"
  )

RNA_counts_annot <- RNA_counts_rep %>%
  left_join(ATL_GO_collapsed, by =c("Geneid" = "Gene")) %>%
  left_join(
    TAIR_Annot %>%
      select(Gene_ID, Gene_name),
    by = c("TAIR_ID" = "Gene_ID")
  )



# GO ANALYSIS ---------

#library dependencies
library(UpSetR)
library(clusterProfiler)
library(pathview)
library(enrichplot)
library(ggplot2)
library(stringr)
library(dplyr)

# --Setting up clusterProfiler and annotations------
# Specify organism annotation database (Arabidopsis thaliana, TAIR)
organism = "org.At.tair.db"



# Load the organism annotation package dynamically
library(organism, character.only = TRUE)

# Display available key types for the annotation database
keytypes(org.At.tair.db)



# D8 vs D1 enrichment analysis########
# Read in DESeq2-identified significantly differentially expressed genes
D1_vs_D8_DEGs <- Regen_DEGs_TAIR_Annot %>%
  dplyr::filter(Comparison == "D1_vs_D8") %>%
  dplyr::distinct(Gene, .keep_all = TRUE)
head(D1_vs_D8_DEGs)


# Extract gene IDs for enrichment analysis
D1_vs_D8_SIG_Genes <- D1_vs_D8_DEGs$TAIR_ID

# Perform GO over-representation analysis (ORA)
# - BP ontology
# - TAIR gene identifiers
# - FDR-adjusted p-value cutoff of 0.05
D1_vs_D8_ORA<- enrichGO(
  gene = D1_vs_D8_SIG_Genes,
  ont = "BP",
  keyType = "TAIR",
  OrgDb = organism,
  pvalueCutoff = 0.05,
  pAdjustMethod = "fdr",
  readable = FALSE
)

# Preview enrichment results
view(D1_vs_D8_ORA)

#Convert enrichResult object to a data frame
D1_vs_D8_GO_Results <- as.data.frame(D1_vs_D8_ORA)


#####Visualize GO enrichment with dotplot------
# Create a dotplot of the top 15 enriched GO BP terms
dotplot <- dotplot(D1_vs_D8_ORA, showCategory = 15) +
  ggtitle("D8vsD1, BP Ontology")

# Center the title and reduce y-axis text size
dotplot +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_text(size = 8)
  )



# Define GO categories of interest for visualization
D1_vs_D8_GO_Top15 <- D1_vs_D8_GO_Results %>%
  dplyr::arrange(desc(FoldEnrichment)) %>%
  dplyr::slice_head(n = 15)
D1_vs_D8_GO_Top15_Descriptions <- D1_vs_D8_GO_Top15$Description


D1_vs_D8_ORA_top15 <- D1_vs_D8_ORA
D1_vs_D8_ORA_top15@result <- D1_vs_D8_ORA@result %>%
  dplyr::filter(Description %in% D1_vs_D8_GO_Top15$Description)

dotplot(D1_vs_D8_ORA_top15 ,showCategory=15) +
  ggtitle("Top 15 GO terms by Fold Enrichment") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_text(size = 8)
  )


# Create a category–gene network plot (cnetplot)
# Nodes represent genes and GO terms
# Edge connections indicate gene membership in GO categories
cnetplot(
  D1_vs_D8_ORA,
  showCategory = D1_vs_D8_GO_Top15_Descriptions,
  layout = igraph::layout_nicely,
  node_label = "category",
  foldChange = D1_vs_D8_ORA@result$FoldEnrichment
)






#Visualization of genes within our Timecourse --------
#### Map DEGs to actual timepoints (convert contrasts → biological timepoints) ----

DEG_map <- Regen_DEGs %>%
  mutate(
    # Translate contrast labels into the timepoint where change occurs
    Timepoint = case_when(
      Comparison == "D1_vs_D8"  ~ "8D",
      Comparison == "D1_vs_D14" ~ "14D",
      Comparison == "D1_vs_D28" ~ "28D"
    )
  ) %>%
  select(Gene, Timepoint) %>%   # keep only what’s needed for mapping
  distinct() %>%               # ensure one entry per Gene–Timepoint pair
  mutate(DEG = TRUE)           # flag used later to mark DEGs in plots



####Prepping Dataframe to Plot Genes in Timecourse ---------

# Extract normalized expression matrix
vsd_mat <- assay(dds_vst)

# Convert FULL matrix → tidy format (no pre-filtering)
Tidy_Regen_df <- as.data.frame(vsd_mat) %>%
  tibble::rownames_to_column("Gene") %>%
  tidyr::pivot_longer(-Gene, names_to = "Sample", values_to = "Expression")

# Attach metadata
Tidy_Regen_Meta <- as.data.frame(colData(dds_vst)) %>%
  tibble::rownames_to_column("Sample")

Tidy_Regen_df <- Tidy_Regen_df %>%
  left_join(Tidy_Regen_Meta, by = "Sample")

# Annotate genes (optional but useful for function search)
Tidy_Regen_df <- Tidy_Regen_df %>%
  left_join(ATL_Annotation, by = "Gene")

Tidy_Regen_df <- Tidy_Regen_df %>%
  left_join(DEG_map, by = c("Gene","Timepoint")) %>%
  mutate(DEG = ifelse(is.na(DEG), FALSE, DEG))

# Enforce timepoint order
Tidy_Regen_df$Timepoint <- factor(Tidy_Regen_df$Timepoint,
                              levels = c("1D","8D","14D","28D"))


# Base Plot Tidy Regen_df --------
Tidy_Regen_df %>%
  filter(Gene == "Soltu.Atl_v3.01_1G000020.2") %>%
  ggplot(aes(x = Timepoint, y = Expression, group = 1)) +
  geom_line() +
  geom_point(size = 3) +
  labs(
    title = "Soltu.Atl_v3.01_1G000020.2",
    x = "Timepoint",
    y = "Expression (VST)"
  )+
  theme_light()

# Plot Function ------

plot_gene_expression <- function(Tidy_Regen_df, Regen_DEGs, query,
                                 direction = c("all", "up", "down")) {
  
  direction <- match.arg(direction)
  
  # ---- FILTER ----
  if (any(query %in% Tidy_Regen_df$Gene)) {
    df <- Tidy_Regen_df %>% dplyr::filter(Gene %in% query)
  } else if (length(query) == 1) {
    df <- Tidy_Regen_df %>%
      dplyr::filter(grepl(query, Functional_Annotation, ignore.case = TRUE))
  } else {
    stop("Query must be either valid gene IDs or a single annotation string")
  }
  
  if (nrow(df) == 0) stop("No matching genes or annotations found.")
  
  # ---- JOIN log2FC ----
  Regen_DEGs <- Regen_DEGs %>%
    dplyr::distinct(Gene, .keep_all = TRUE)
  
  df <- df %>%
    dplyr::left_join(
      Regen_DEGs %>% dplyr::select(Gene, log2FoldChange),
      by = "Gene"
    )
  

  
  # ---- Direction ----
  df <- df %>%
    dplyr::mutate(
      DEGs = dplyr::case_when(
        DEG & !is.na(log2FoldChange) & log2FoldChange > 0 ~ "Up",
        DEG & !is.na(log2FoldChange) & log2FoldChange < 0 ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  plot_data <-df %>%
    dplyr::select(Gene, Functional_Annotation, Timepoint, Expression, DEG)

  # ---- DYNAMIC LEGEND LOGIC ----
  
  n_genes <- length(unique(df$Gene))
  
  if (n_genes <= 10) {
    
    #  NO filtering — keep all genes
    df <- df %>%
      dplyr::mutate(
        Gene_Label = paste0(
          Gene, " - ",
          stringr::str_trunc(Functional_Annotation, 40)
        )
      )
    
    legend_breaks <- unique(df$Gene_Label)
    
  } else {
    
    #  ONLY apply this when >10 genes
    top_genes <- df %>%
      dplyr::filter(DEG) %>%
      dplyr::group_by(Gene) %>%
      dplyr::summarise(
        max_fc = ifelse(
          all(is.na(log2FoldChange)),
          NA_real_,
          max(abs(log2FoldChange), na.rm = TRUE)
        ),
        .groups = "drop"
      ) %>%
      dplyr::arrange(desc(max_fc)) %>%
      dplyr::slice_head(n = 10)
    
    df <- df %>%
      dplyr::mutate(
        Gene_Label = ifelse(
          Gene %in% top_genes$Gene,
          paste0(
            Gene, " - ",
            stringr::str_trunc(Functional_Annotation, 40)
          ),
          "Other"
        )
      )
    
    legend_breaks <- paste0(
      top_genes$Gene, " - ",
      stringr::str_trunc(
        df$Functional_Annotation[match(top_genes$Gene, df$Gene)],
        40
      )
    )
  } 
  # ---- Hover ----
  df$hover <- paste(
    "Gene:", df$Gene,
    "<br>Function:", df$Functional_Annotation,
    "<br>Timepoint:", df$Timepoint,
    "<br>Expression:", round(df$Expression, 2),
    "<br>log2FC:", round(df$log2FoldChange, 2)
  )
  
  # ---- DEG subset ----
  df_deg <- df %>%
    dplyr::group_by(Gene) %>%
    dplyr::filter(any(DEG)) %>%
    dplyr::ungroup()
  
  # ---- Direction filter ----
  if (direction == "up") {
    genes_keep <- df %>%
      dplyr::group_by(Gene) %>%
      dplyr::filter(any(DEGs == "Up")) %>%
      dplyr::pull(Gene) %>%
      unique()
    
    df_deg <- df %>% dplyr::filter(Gene %in% genes_keep)
    
  } else if (direction == "down") {
    genes_keep <- df %>%
      dplyr::group_by(Gene) %>%
      dplyr::filter(any(DEGs == "Down")) %>%
      dplyr::pull(Gene) %>%
      unique()
    
    df_deg <- df %>% dplyr::filter(Gene %in% genes_keep)
  }
  
  deg_plot_data <- df_deg %>%
    dplyr::select(Gene, Functional_Annotation, Timepoint, Expression, DEG)
  
  # ---- PLOTS ----
  
  p1 <- ggplot(df, aes(Timepoint, Expression,
                       group = Gene,
                       color = Functional_Annotation)) +
    geom_line() +
    geom_point(aes(fill = DEGs),
               shape = 21, size = 3) +
    scale_fill_manual(values = c("Up"="green","Down"="red","NS"="white")) +
    theme_light()
  
  # ---- P2 (FIXED LEGEND) ----
  p2 <- ggplot(df, aes(Timepoint, Expression,
                       group = Gene,
                       color = Gene_Label)) +
    geom_line() +
    geom_point(aes(fill = DEGs),
               shape = 21, size = 3) +
    scale_fill_manual(values = c("Up"="green","Down"="red","NS"="white")) +
    scale_color_discrete(breaks = legend_breaks)+
    theme_light() +
    labs(title = paste("Gene view:", paste(query, collapse = ", ")))
  
  # ---- P3 (FIXED LEGEND) ----
  if (nrow(df_deg) == 0) {
    p3 <- NULL
    p3_int <- NULL
  } else {
    
    p3 <- ggplot(df_deg, aes(Timepoint, Expression,
                             group = Gene,
                             color = Gene_Label)) +
      geom_line() +
      geom_point(aes(fill = DEGs),
                 shape = 21, size = 3) +
      scale_fill_manual(values = c("Up"="green","Down"="red","NS"="white")) +
      scale_color_discrete(breaks = legend_breaks)+
      theme_light() +
      labs(title = paste("DEG-only genes"))
    
    p3_int <- ggplotly(p3 + aes(text = hover), tooltip = "text")
  }
  
  # ---- INTERACTIVE ----
  p1_int <- ggplotly(p1 + aes(text = hover), tooltip = "text")
  p2_int <- ggplotly(p2 + aes(text = hover), tooltip = "text")
  
  if (!is.null(p3)) {
    p3_int <- ggplotly(p3 + aes(text = hover), tooltip = "text")
  }
  
  return(list(
    static_function_plot = p1,
    static_gene_plot = p2,
    interactive_function_plot = p1_int,
    interactive_gene_plot = p2_int,
    deg_only_plot = p3,
    deg_only_interactive_plot = p3_int,
    plot_data =plot_data,
    deg_plot_data = deg_plot_data 
  ))
}


#Plot Example 1: Using Annotation Functionality -------
TCP_plot <- plot_gene_expression(Tidy_Regen_df,Regen_DEGs, "TCP")
TCP_plot$static_function_plot
#TCP_plot$interactive_function_plot
TCP_plot$static_gene_plot
#TCP_plot$interactive_gene_plot
TCP_plot$deg_only_plot
TCP_plot$interactive_gene_plot

TCP_plot_df<-as.data.frame(TCP_plot$plot_data)
TCP_plot_DEGs_df <- as.data.frame(TCP_plot$deg_plot_data)

TCP_plot_up<- plot_gene_expression(Tidy_Regen_df,Regen_DEGs, "TCP",direction= "up")
TCP_plot_up$static_function_plot
#TCP_plot_up$interactive_function_plot
TCP_plot_up$static_gene_plot
#TCP_plot_up$interactive_gene_plot
TCP_plot_up$deg_only_plot
#TCP_plot_up$interactive_gene_plot

TCP_plot_down<- plot_gene_expression(Tidy_Regen_df,Regen_DEGs, "TCP",direction= "down")
TCP_plot_down$static_function_plot
#TCP_plot_down$interactive_function_plot
TCP_plot_down$static_gene_plot
#TCP_plot_down$interactive_gene_plot
TCP_plot_down$deg_only_plot
#TCP_plot_down$interactive_gene_plot


# Plot Example 2: Plot specific genes of interest (manual gene list) -------

# Define curated gene list (e.g., tuber-related candidates)
Tuber_genes <- c(
  "Soltu.Atl_v3.05_1G024440.1",
  "Soltu.Atl_v3.05_2G024320.1",
  "Soltu.Atl_v3.05_2G024850.1",
  "Soltu.Atl_v3.05_1G022420.1",
  "Soltu.Atl_v3.05_1G022430.1",
  "Soltu.Atl_v3.05_1G022450.1",
  "Soltu.Atl_v3.05_2G022000.1",
  "Soltu.Atl_v3.05_2G022020.1",
  "Soltu.Atl_v3.05_3G022530.1",
  "Soltu.Atl_v3.05_3G022510.1",
  "Soltu.Atl_v3.05_4G021700.1",
  "Soltu.Atl_v3.05_4G021720.1",
  "Soltu.Atl_v3.06_2G021580.1",
  "Soltu.Atl_v3.06_3G015520.1",
  "Soltu.Atl_v3.06_1G010970.2",
  "Soltu.Atl_v3.06_2G022190.2",
  "Soltu.Atl_v3.06_3G016130.1",
  "Soltu.Atl_v3.06_4G015370.2",
  "Soltu.Atl_v3.06_1G014050.1",
  "Soltu.Atl_v3.06_2G026030.4",
  "Soltu.Atl_v3.06_3G019920.3",
  "Soltu.Atl_v3.02_1G027420.1",
  "Soltu.Atl_v3.02_2G029850.1",
  "Soltu.Atl_v3.02_4G022030.1",
  "Soltu.Atl_v3.02_1G027430.1",
  "Soltu.Atl_v3.02_2G029870.1",
  "Soltu.Atl_v3.02_2G030560.1",
  "Soltu.Atl_v3.02_2G030580.1",
  "Soltu.Atl_v3.02_4G022040.1"
)

# Ensure proper format (character vector, remove duplicates)
Tuber_genes <- as.character(Tuber_genes)
Tuber_genes <- unique(Tuber_genes)


# Plot full time-course expression for selected genes
Tuber_genes_plots <- plot_gene_expression(Tidy_Regen_df, Regen_DEGs,Tuber_genes)

# Static plots (clean for publication/export)
Tuber_genes_plots$static_function_plot
Tuber_genes_plots$static_gene_plot
Tuber_genes_plots$deg_only_plot

# Interactive plots (exploratory analysis / hover details)
Tuber_genes_plots$interactive_function_plot
Tuber_genes_plots$interactive_gene_plot
Tuber_genes_plots$deg_only_interactive_plot



# Plot Example 3: Use external annotation table (TF families) -------

# Load transcription factor annotation table
Potato_TF_geneids<-read.csv("GHC_STuberosum_TF_annotation.csv",sep=',')

# Keep only relevant columns and standardize naming
Potato_TF_geneids <- Potato_TF_geneids[, c("Subject_accession", "Family")]
colnames(Potato_TF_geneids)<-c("Gene", "Family")

# Extract all TF gene IDs
Potato_TF<-c(Potato_TF_geneids$Gene)

# Plot ALL transcription factors (broad overview)
TF_plots <- plot_gene_expression(Tidy_Regen_df, Regen_DEGs,Potato_TF)
TF_plots$static_function_plot

# Export DEG subset used in plot for downstream analysis
write.table(degs_plot_df,"TF_DEGs_Regen_Trial_DF.csv",sep="\t")


# Focus on a specific TF family (e.g., LBD) ----
TF_Family_name <- "LBD"   

# Subset genes belonging to selected TF family
TF_Family_genes <- Potato_TF_geneids %>%
  dplyr::filter(Family == TF_Family_name) %>%
  dplyr::pull(Gene)

# Keep only genes present in expression dataset
TF_Family_genes <- intersect(TF_Family_genes, Tidy_Regen_df$Gene)

# Plot expression for selected TF family (more interpretable than plotting all TFs at once)
TF_plots <- plot_gene_expression( Tidy_Regen_df, Regen_DEGs, TF_Family_genes)

#Gene static visualization 
TF_plots$static_function_plot
TF_plots$static_gene_plot

# DEG-focused visualization (only genes with differential signal)
TF_plots$deg_only_plot
TF_plots$deg_only_interactive_plot

# Extract DEG plot data for export or further analysis
TF_plots_df <- as.data.frame(TF_plots$deg_plot_data)


# Plot only downregulated TFs (direction-specific filtering)
TF_plots <- plot_gene_expression(Tidy_Regen_df, Regen_DEGs,Potato_TF,direction = "down")

# DEG-focused visualization (only genes with differential signal)
TF_plots$deg_only_plot

# Plot only downregulated TFs (direction-specific filtering)
TF_plots <- plot_gene_expression(Tidy_Regen_df, Regen_DEGs,Potato_TF,direction = "up")

# DEG-focused visualization (only genes with differential signal)
TF_plots$deg_only_plot

