#installing edgeR library
BiocManager::install("edgeR")
library(edgeR)#importing library

#importing the data, making Geneid the first col, making header true
counts = read.table("counts.txt",row.names = 1, header=T)
counts<-counts[, -(1:5)]# removing irrelevant columns
dge <- DGEList(counts = counts)

# Convert to CPM
cpm_counts <- cpm(dge)

# Keep genes with at least 1 CPM for all samples
keep <- rowSums(cpm_counts >= 1) == 3
dge <- dge[keep, , keep.lib.sizes=FALSE]

# TMM normalisation
dge <- calcNormFactors(dge, method="TMM")
###############################################
# Define columns for each group (updating sample names from metadata)
# Rename columns of counts
colnames(counts) <- c("M_G1", "S", "G2")
colnames(counts)

norm_counts <- cpm(dge, normalized.lib.sizes = TRUE)

# Mean log2 expression per group
mean_MG1 <- rowMeans(log2(norm_counts[, M_G1,drop = FALSE] + 1))
mean_S   <- rowMeans(log2(norm_counts[, S,drop = FALSE] + 1))
mean_G2  <- rowMeans(log2(norm_counts[, G2,drop = FALSE] + 1))

# Calculate log2 fold changes
log2FC_MG1_vs_S  <- mean_MG1 - mean_S
log2FC_S_vs_G2   <- mean_S   - mean_G2
log2FC_G2_vs_MG1 <- mean_G2  - mean_MG1

# Combine into one results table, keeping the same index as norm_counts
results <- data.frame(
  log2FC_MG1_vs_S  = log2FC_MG1_vs_S,
  log2FC_S_vs_G2   = log2FC_S_vs_G2,
  log2FC_G2_vs_MG1 = log2FC_G2_vs_MG1
)

# Preserve the gene IDs as rownames
rownames(results) <- rownames(norm_counts)

head(results)
##############
library(limma)

# --- Define sample groups ---
# Replace these with the actual sample labels in your counts table
group <- factor(c(rep("M_G1", 1), rep( "S", 1), rep("G2", 1)))

# --- Design matrix ---
design<- model.matrix(~0 + group)
colnames(design) <- levels(group)

# --- Transform counts with voom ---
v <- voom(dge, design, plot = FALSE)

# --- Fit linear model ---
fit <- lmFit(v, design)

# --- Define contrasts ---
contrasts <- makeContrasts(
  M_G1_vs_S  = M_G1 - S,
  S_vs_G2    = S - G2,
  G2_vs_M_G1 = G2 - M_G1,
  levels = design
)

fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2) # No residual degrees of freedom in linear model fits
