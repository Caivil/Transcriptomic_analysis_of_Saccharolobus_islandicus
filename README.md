# Transcriptomic_analysis_of_*Saccharolobus_islandicus*
## 1.Data source
The experimental design of the study focused on understanding the gene expression of Saccharolobus islandicus relative to its different cell-cycle stages.  Cultures of S. islandicus were made to begin the cell cycle at once by a 6-hour treatment with acetic acid. After the removal of the acetic acid, the cells resumed their cell cycle. Total RNA was extracted from 3 samples at three specific time points after synchronization:
- sample 1:  2 hours and 30 minutes (in M-G1 phase)
- sample 2: 4 hours (in S phase)
- sample 3: 6 hours( in G2 phase)
## 2.Obtaining the data
The reads (Sequenced with Illumina NextSeq 2000) were obtaiend from the NCBI database along with a ```.fasta``` reference and ```gunzip *.fastq.gz``` was used to unzip the file. 
## 3.Workflow
This workflow explains the sequence processing of subsamples
### In Bash
#### 3.1 Mapping reads to reference 
```module load bowtie2-2.4.1, module load samtools-1.18,module load subread-2.0.8```  was ran to load the tool into the enivironmnet.
###### Loop through each sample
```samples=(2h30_13_S26_R1_001.fastq  4h_1_S31_R1_001.fastq  6h_10_S53_R1_001.fastq)
for i in "${samples[@]}"; do
    bowtie2 -x ref_index -U "${i}.fastq" -S "${i}.sam"
    samtools view -b "${i}.sam" > "${i}.bam"
    rm "${i}.sam"
done
```
#### 3.2 FeatureCounts for quantifying reads
```
featureCounts -a sequence.gtf -F GTF -o counts_2.txt -T 14 *.bam
```
### In R
The ```edgeR library``` was instealed and imported
```
Installing edgeR library
BiocManager::install("edgeR")
library(edgeR)#importing library
```
### 3.3 importing the data, making Geneid the first col, making header true and removing irrelevant columns
```
counts = read.table("counts.txt",row.names = 1, header=T)
counts<-counts[, -(1:5)]
dge <- DGEList(counts = counts)
```
### 3.4 Convert to CPM
```
cpm_counts <- cpm(dge)
```
### 3.5 Keep genes with at least 1 CPM for all samples
```
keep <- rowSums(cpm_counts >= 1) == 3
dge <- dge[keep, , keep.lib.sizes=FALSE]
```
### 3.6 TMM normalisation, define columns for each group (updating sample names from metadata) and rename columns of counts
```
dge <- calcNormFactors(dge, method="TMM")
colnames(counts) <- c("M_G1", "S", "G2")
colnames(counts)
norm_counts <- cpm(dge, normalized.lib.sizes = TRUE)
```
### 3.7 Mean log2 expression per group was determined
```
mean_MG1 <- rowMeans(log2(norm_counts[, M_G1,drop = FALSE] + 1))
mean_S   <- rowMeans(log2(norm_counts[, S,drop = FALSE] + 1))
mean_G2  <- rowMeans(log2(norm_counts[, G2,drop = FALSE] + 1))
```
### 3.8 Calculate log2 fold changes
```
log2FC_MG1_vs_S  <- mean_MG1 - mean_S
log2FC_S_vs_G2   <- mean_S   - mean_G2
log2FC_G2_vs_MG1 <- mean_G2  - mean_MG1
```
### 3.9 Combine into one results table, keeping the same index as norm_counts
```
results <- data.frame(
  log2FC_MG1_vs_S  = log2FC_MG1_vs_S,
  log2FC_S_vs_G2   = log2FC_S_vs_G2,
  log2FC_G2_vs_MG1 = log2FC_G2_vs_MG1
)
```
### 3.10 Preserve the gene IDs as rownames
```
rownames(results) <- rownames(norm_counts)
head(results)
```

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



