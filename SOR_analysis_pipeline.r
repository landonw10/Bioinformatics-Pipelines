# data available at: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE223066


# Aggregate raw gene expression counts across all barcodes for each mouse, 
# producing a genes × samples pseudobulk matrix

library(rhdf5)

# Initialize a list to store a raw count vector
all_counts <- list()

# Iterate over each sample
for (i in 1:14) {
  file_path <- paste0("/Visium/spatial ", i, "/filtered_feature_bc_matrix.h5")
  print(paste("Processing sample", i))
  
  # Read data from HDF5 files
  features <- h5read(file_path, "/matrix/features/name")
  data <- h5read(file_path, "/matrix/data")
  indices <- h5read(file_path, "/matrix/indices")  # Indices are positions of the genes in 'features'
  
  # Create an empty vector to hold the aggregated expression for each gene
  count_vector <- rep(0, length(features))
  
  # Aggregate counts for each gene across all locations in the current sample
  for (j in 1:length(data)) {
    gene_index <- indices[j] + 1  # Convert to 1-indexing (h5 file uses 0-indexing)
    count_vector[gene_index] <- count_vector[gene_index] + data[j]
  }
  
  # Add sample(i) matrix of gene counts to list
  count_matrix <- matrix(count_vector, nrow = length(features), ncol = 1, dimnames = list(features, paste0("Sample_", i)))
  all_counts[[i]] <- count_matrix
}

# Combine the aggregated counts into a single matrix (genes x samples)
counts_matrix <- do.call(cbind, all_counts)


# Perform differential expression analysis using edgeR on pseudobulk data,
# comparing HC vs SOR groups based on aggregated gene expression counts.

library(edgeR)

# Define group labels (0 = HC 1 = SOR)
group_labels <- rep(0, 14)
group_labels[c(1, 2, 3, 4, 10, 12, 14)] <- 1

# Create a DGEList object for edgeR
dge <- DGEList(counts = counts_matrix)

# Filter out genes with fewer than 10 counts in any sample
keep <- rowSums(dge$counts >= 10) >= 1 
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Normalize the data
dge <- calcNormFactors(dge)

# Assign groups from group_labels
dge$samples$group <- group_labels

# Design matrix for differential expression analysis
design <- model.matrix(~ group, data = dge$samples)

# Estimate dispersion
dge <- estimateDisp(dge, design)

# Fit the model
fit <- glmFit(dge, design)

# Perform likelihood ratio test
lrt <- glmLRT(fit)

# Get the top differential expressed genes
top_genes <- topTags(lrt, n = Inf)

# Report the number of up- and down-regulated genes
upregulated <- sum(top_genes$table$logFC > 0.5 & top_genes$table$FDR < 0.05)
downregulated <- sum(top_genes$table$logFC < -0.5 & top_genes$table$FDR < 0.05)

cat("Number of up-regulated genes:", upregulated, "\n")
cat("Number of down-regulated genes:", downregulated, "\n")

# Report the top 5 most up- and down-regulated genes
top_upregulated <- head(top_genes$table[top_genes$table$logFC > 0.5 & top_genes$table$FDR < 0.05, ], 5)
top_downregulated <- head(top_genes$table[top_genes$table$logFC < -0.5 & top_genes$table$FDR < 0.05, ], 5)

cat("Top 5 up-regulated genes:\n")
print(rownames(top_upregulated))

cat("Top 5 down-regulated genes:\n")
print(rownames(top_downregulated))


# Perform GO enrichment analysis for molecular function (MF) using upregulated and downregulated genes

library(clusterProfiler)
library(org.Mm.eg.db)

# Convert upregulated and downregulated gene lists to Ensembl IDs
up <- c(rownames(top_upregulated))
down <- c(rownames(top_downregulated))
up_entrez <- bitr(up, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)$ENTREZID
down_entrez <- bitr(down, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)$ENTREZID

# Perform GO enrichment analysis for molecular function
up_go <- enrichGO(gene = up_entrez, OrgDb = org.Mm.eg.db, ont = "MF", pvalueCutoff = 0.05)
down_go <- enrichGO(gene = down_entrez, OrgDb = org.Mm.eg.db, ont = "MF", pvalueCutoff = 0.05)
top_up_go <- head(up_go, 5)
top_down_go <- head(down_go, 5)

# Report results with biological descriptions
cat("Top 5 GO molecular function terms enriched in upregulated genes:\n")
for (i in seq_len(nrow(top_up_go))) {
  cat(paste(top_up_go$Description[i], "\n"))
}

cat("\nTop 5 GO molecular function terms enriched in downregulated genes:\n")
for (i in seq_len(nrow(top_down_go))) {
  cat(paste(top_down_go$Description[i], "\n"))
}


# Perform SOR vs HC differential expression analysis on count bulk RNA-seq data from tximport, 
# using edgeR's exactTest method and retrieving the top 5 differentially expressed genes.

library(biomaRt)

file <- "/GSE223066_tximport-counts_learning1.txt"
data <- read.table(file, header = TRUE, row.names = 1)
data_rounded <- round(data)

# Define the group labels: "HC" as control and "SOR" as experimental
group_labels <- c(rep("HC", 4), rep("SOR", 4))

# Create a DGEList object
dge_data <- DGEList(counts = data_rounded)
dge_data$samples$group <- group_labels

# Normalize the data (library size normalization)
dge_data <- calcNormFactors(dge_data)

# Perform the differential expression analysis using the exactTest function
dge_data <- estimateDisp(dge_data) 
result <- exactTest(dge_data)

# Extract the top 5 genes with the smallest p-values as Ensembl IDs
top_genes <- topTags(result, n = 5)
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
top_genes_ensg <- rownames(top_genes$table)
top_genes_symbols <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                           filters = "ensembl_gene_id", values = top_genes_ensg, mart = ensembl)
top_genes_with_symbols <- merge(top_genes$table, top_genes_symbols, by.x = "row.names", by.y = "ensembl_gene_id")

# Print the gene symbols of the top genes
print(top_genes)
print(top_genes_with_symbols$external_gene_name)


# Compare log-fold changes between pseudobulk and bulk RNA-seq datasets to determine correlation

# Get gene symbols from the pseudobulk dataset as Ensembl IDs
pseudobulk_genes <- rownames(dge$counts)
pseudobulk_ensembl <- bitr(pseudobulk_genes, fromType = "SYMBOL", toType = "ENSEMBL", OrgDb = org.Mm.eg.db)

# Identify gene symbols that map to multiple Ensembl IDs
duplicated_symbols <- pseudobulk_ensembl$SYMBOL[duplicated(pseudobulk_ensembl$SYMBOL) | 
                                                  duplicated(pseudobulk_ensembl$SYMBOL, fromLast = TRUE)]

# Remove duplicates and double matches
pseudobulk_ensembl_unique <- pseudobulk_ensembl[!pseudobulk_ensembl$SYMBOL %in% duplicated_symbols, ]
pseudobulk_ensembl_unique <- pseudobulk_ensembl_unique[!duplicated(pseudobulk_ensembl_unique$ENSEMBL), ]

# Use Ensembl IDs from bulk data directly
bulk_ensembl <- rownames(data_rounded)

# Find common Ensembl IDs between pseudobulk and bulk datasets
common_ensembl <- intersect(pseudobulk_ensembl_unique$ENSEMBL, bulk_ensembl)
num_common_genes <- length(common_ensembl)
cat("Number of symbols that appear in both experiments:", num_common_genes, "\n")

# Map pseudobulk gene symbols to Ensembl IDs
pseudobulk_lfc <- lrt$table[match(pseudobulk_ensembl_unique$SYMBOL[pseudobulk_ensembl_unique$ENSEMBL %in% common_ensembl], rownames(lrt$table)), "logFC"]

# Extract log-fold changes for common genes from bulk data
bulk_lfc <- result$table[common_ensembl, "logFC"]

# Perform a correlation analysis
correlation <- cor(pseudobulk_lfc, bulk_lfc, method = "pearson")

# Perform a linear regression analysis
regression <- lm(bulk_lfc ~ pseudobulk_lfc)

# Report the correlation coefficient and regression summary
cat("Correlation coefficient between pseudobulk and bulk log-fold changes:", correlation, "\n")
summary(regression)

# Plot a scatter plot of the log-fold changes for common genes
plot(pseudobulk_lfc, bulk_lfc, 
     xlab = "Pseudobulk Log-Fold Change", 
     ylab = "Bulk Log-Fold Change", 
     main = "Log-Fold Change Comparison: Pseudobulk vs. Bulk", 
     pch = 1,
     col = rgb(0, 0, 1, 0.5))  # Blue with some transparency
abline(lm(bulk_lfc ~ pseudobulk_lfc), col = "red", lwd = 2)
grid()