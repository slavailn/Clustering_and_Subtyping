# ConsensusClusterPlus (CCP) Workflow
A short tutorial dedicated to identification of tumor molecular subtypes based on gene expression profiles.   
---
## Prerequisites

> Installs and loads all R/Bioconductor packages used later in the script.

```r
# Install once (skip if already installed)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "ConsensusClusterPlus",     # the star of the show
  "ComplexHeatmap", "circlize", # prettier consensus matrices
  "cluster", "factoextra"       # silhouette metrics & plots
))

install.packages(c("tidyverse", "RColorBrewer", "reshape2", "vcd"))

# Load libraries for this session
library(tidyverse)
library(ConsensusClusterPlus)
library(ComplexHeatmap)
library(circlize)
library(cluster)
library(factoextra)
```
---

## 1. Load input data

> Read the expression matrix (genes × samples) and the sample‑level metadata, then forces the two to share the same sample order.
> We assume that the expression data was normalized and variance stabilized.
```r
expr <- read.csv("expression_matrix_norm_vsd.csv", header = T)    # normalized and variance stabilized data
metadata  <- read.csv("sample_metadata.csv", header = T)  # load metadata 

expr <- column_to_rownames(expr, var = "Gene")
expr  <- expr[ , metadata$SampleID ] # enforce order match
```

---

## 2.  Feature selection
Clustering is performed based on a certain number of most variable genes. 
Typically, researchers select anywhere between 500 - 5000 top most variable genes to include in the analysis.
Published bulk-tumour studies and the ConsensusClusterPlus vignette pick a value in the ~500 – 5 000 most-variable genes band and tune 
it following some stability metrics.  

Gene subset	When it works well	Typical cohort size
500 – 1 000:	Small cohorts (< 50 samples) or noisy data.
1 500 – 3 000	Mid-sized studies (50-200 samples).	Most TCGA-scale papers land here.
≈ 5 000	Large, well-normalized datasets (> 200 samples) or legacy microarrays	You have power to support more features.

Rule of thumb
Aim for “features ≤ 20 × samples”. With 100 tumours, 2 000 genes is safely inside that bound.
Best practice (and what the workflow does)
    Start at 500 genes.
    Run CCP, check PAC + CDF elbow.
    Try 500, 1000, 2 000, 3 000, and 5000 genes.

    Keep the smallest N that still gives PAC ≤ 0.15 and a clear delta area plateau.

```r
# Select top 1000 genes by median absolute deviation (MAD)
disp <- apply(exp, 1, mad)
ranked <- exp[order(disp, decreasing = T),]
head(ranked)
expr1000 <- as.matrix(ranked[1:1000,])

```

---

## 3. Run ConsensusClusterPlus (hierarchical + Pearson)

> Bootstraps the dataset 1 000 times, each time clustering 80 % of samples with hierarchical clustering on 1 000 genes, 
using 1 – Pearson correlation** as distance. The distance measure was "pearson" and clustering algorithm was "hc".
> Try CCP with the following cutoffs of top most variable genes: 500, 500, 1000, 2000, 3000, 5000

```r
set.seed(42)   # reproducible
# Suggested cutoffs
cutoffs <- c(500, 1000, 2000, 3000, 5000)

# For example run CCP with top 500 genes.
# Try all of the cutoffs suggested above, save the results of each run
cc <- 
  ConsensusClusterPlus(as.matrix(ranked[1:500,]), 
                        maxK = 15, reps = 1000,
                        pItem = 0.8, pFeature = 1,
                        clusterAlg = "hc", distance = "pearson",
                        seed = 1, plot = "png", 
                        title = "CC_plots_top500")
save(cc, file = "CC_top500.RData")

```

---

## 4. Pick the optimal **N** top genes that create the most stable cluster structure. 
> Calculates PAC (Proportion of Ambiguously Clustered pairs) for each number of clusters K; chooses the smallest K where PAC ≤ 0.15 and the delta‑area curve flattens.
> Balances cluster stability (low PAC) with parsimony (few clusters).

```r
# Calculate PAC
lower <- 0.10   # lower consensus threshold
upper <- 0.90   # upper consensus threshold

# cc is the list returned by ConsensusClusterPlus
# it is indexed by K = 2,3,...,maxK (position 1 is K=2)
# pac_one is a function helper that converts a consensus matrix to PAC
pac_one <- function(consmat, lo = lower, hi = upper) {
  # take the upper triangle without the diagonal
  v <- consmat[upper.tri(consmat, diag = FALSE)]
  mean(v > lo & v < hi)   # proportion of ambiguous pairs
}
Ks  <- 2:length(cc) # the K values actually present
PAC <- sapply(Ks, function(k) pac_one(cc[[k]]$consensusMatrix))

# put results in a tidy data.frame
pac_df <- data.frame(K = Ks, PAC = PAC)
pac_df

# Top 6 rows of PAC data frame
# K       PAC
# 2 0.0155642
# 3 0.4712997
# 4 0.1437276
# 5 0.1639067
# 6 0.1783247

# Interpret PAC table
# PAC (Proportion of Ambiguously Clustered pairs) measures the fraction of sample pairs whose 
# consensus values fall inside an “uncertain” range (typically 0.1–0.9).
# 0 in every pair is either always together or never together (perfectly stable).
# Values < 0.05 are considered excellent, 0.05–0.10 good, 0.10–0.15 acceptable, and > 0.15 increasingly unstable.
# K = 2 is the only highly stable solution (PAC ≈ 0.0155642).
# Jumping to K = 3 PAC increases to 0.471, which indicates very unstable cluster structure.
# The lowest PAC is 0.143 observed at K = 4.
# Practical guidance
# If biological or clinical considerations demand ≥ 3 clusters, you should:
# Re-evaluate feature selection (e.g., use more variable genes or a different distance metric) to see if PAC for K = 3 can be pushed below ~0.15.
# Otherwise, consider alternative strategies such as accepting the very stable 2-cluster split and performing sub-clustering within each of those two groups.
# If interpretability is not compromised, K = 2 is statistically the most defensible choice for this feature set.

# In this case K = 2 produces lowest PAC, but this split is not biologically relevant. On the other hand, K = 4 produces acceptable 
# cluster stablity and is in good agreement with histological types of tumors in the dataset tested.  

```
**Note!** After examining PAC tables generated with top 500, 1000, 2000, 3000, and 5000 most variable genes, the lowest PAC outside of K > 2 
was observed in K = 4 with top 1000 genes. Further fine-tuning was done with 1000 most variable genes (based on MAD). I rejected K = 2 cluster structure based on 
biological irrelevance, since dataset contained multiple tumor types.

## 5. Test various combinations of clustering parameters
Create parameters grid to form different combinations of clustering algorithms and distance metrics.
```r
# Create parameters grid
algs   <- c("hc", "hc", "hc", "km", "pam")
dists  <- c("pearson","spearman","euclidean" ,"euclidean","pearson")
labels <- LETTERS[1:6]

results <- vector("list", length(algs))
names(results) <- labels

# Iterate over combinations of parameters and apply consensus clustering
# We use top 1000 variable genes since we established that at this gene number
# threshold provides most stable cluster structure
for (i in seq_along(algs)) {
  cat("Running", labels[i], algs[i], dists[i], "...\n")
  mat <- if (dists[i] == "euclidean") t(scale(t(expr1000))) else expr1000
  results[[i]] <- ConsensusClusterPlus(as.matrix(mat),
                                       maxK = 7, reps = 1000,
                                       pItem = 0.80, pFeature = 1,
                                       clusterAlg = algs[i], distance = dists[i],
                                       seed = 42, plot = NULL)
}

# Calculate PAC
PAC <- function(M, lo = 0.10, hi = 0.90) {
  F <- ecdf(M[upper.tri(M)]); F(hi) - F(lo)
}

# Aggregate PAC scores for each set of parameter combinations
score_tbl <- do.call(rbind, lapply(labels, \(lab) {
  pac <- sapply(results[[lab]][2:5], \(x) PAC(x$consensusMatrix))
  data.frame(Method = lab, K = 2:5, PAC = pac)
}))
print(score_tbl)

score_tbl$method_combo <- paste(algs, dists, sep = "_")
write.csv(score_tbl, file = "PAC_grid_test.csv", row.names = F)

# Comnbination of distance measure based on pearson correlation and hierarchical clustering still provides the
# best PAC, which remains 0.14 at K = 4

```
**Note!** Searching parameter combinations established 4 cluster structure with "pearson" correlation and hierarchical
clustering as most stable. We will continue the analysis using this structure.

## 6. Extract cluster labels

>  Pulls the final cluster assignment vector for K = 4 and saves it.
>  Downstream analyses (silhouette, histology, survival) need these labels.

```r
clusters <- cc_final[[4]]$consensusClass   # named by SampleID
write_tsv(tibble(SampleID = names(clusters),
                 Cluster   = clusters),
          "Subtypes_K4.tsv")
```

---

## 7. Stability diagnostics

### 6.1 Silhouette

> Measures how well each sample fits within its cluster.
> PAC summarises pairwise stability; silhouette gives a geometry‑based view.

```r
diss_spear <- as.dist(1 - cor(expr_sel, method = "spearman"))

sil_sp <- silhouette(clusters, diss_spear)   # mean ≈ 0.24 (borderline‑good)
```

### 6.2 Plot silhouette (Spearman)

```r
png("silhouette_K4_spearman.png", 1600, 900, res = 180)
fviz_silhouette(sil_sp,
                palette = brewer.pal(4, "Set2"),
                label   = FALSE,
                print.summary = TRUE) +
  labs(title = "Silhouette (K = 4, Spearman distance)")
dev.off()
```
## 8. Integrate clusters with metadata

> Added cluster labels to metadata and tested independence between clusters and categorical variables (histology, stage).
> Validates biological relevance and reveals enrichments.

```r
# Added cluster memberships to the metadata tables
metadata_aug <- metadata %>%
  left_join(read_csv("histology_lookup.csv"), by = "Histology") %>%
  mutate(Cluster = clusters[ SampleID ])

# Chi square test: Cluster × Histology
hist_tab <- table(metadata_aug$Cluster, metadata_aug$Histology_abbr)
hist_chi <- chisq.test(hist_tab)
cramersV <- sqrt(hist_chi$statistic /
                 (sum(hist_tab) * (min(dim(hist_tab)) - 1)))
```

_(Cramer’s V ≈ 0.74 indicates strong association.). Very high significance in Chi square test_

### Plot 100 % stacked bars

```r
ggplot(metadata_aug, aes(factor(Cluster), fill = Histology_abbr)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Cluster (K = 4)", y = "% of samples",
       title = "Histology distribution across clusters") +
  theme_bw()
```
### Identify which histology-cluster pair drives the signal
```r
std_res <- chisq_res$stdres   # standardised Pearson residuals
round(std_res, 2)
```
Use the same type analysis to examine the association between stage and clusters.

---

## 8. Create consensus matrix

> Visualise the 4×4 consensus matrix with annotations (cluster, histology, stage).
> Easy‑to‑read figure for manuscripts; highlights crisp diagonal blocks.

```r
M <- cc_final[[4]]$consensusMatrix
rownames(M) <- colnames(M) <- names(clusters)

ha <- HeatmapAnnotation(
  Cluster   = factor(clusters),
  Histology = metadata_aug$Histology_abbr,
  Stage     = metadata_aug$Stage,
  col = list(
    Cluster   = brewer.pal(4

```



