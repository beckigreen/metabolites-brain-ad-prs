#=================================#
#  0.Metabolic correlates paper:  #
#  relevant analyses (WGCNA/ORA)  #    
#=================================#

## Set up ----
filepath <- ""
setwd(filepath)

### Load libraries ----
library(dplyr)
library(openxlsx)
library(tidyr)
library(WGCNA)

### Load data ----
# Cognition & KNN imputed metabolite data
cog_met <- readRDS("dat/other/cog_met.rds")

# Extract metabolite colnames
metab_names <- cog_met %>%
  select(contains("metabolon_")) %>%
  colnames()

row.names(cog_met) <- cog_met$NSHD_ID

dim(cog_met) 
length(metab_names)

# Select cols required for WGCNA
wgcna_dat <- cog_met %>%
  select(sex, 
         FASTING, 
         BLCLIN09, 
         AGEN09, 
         all_of(metab_names)) 

dim(wgcna_dat)
row.names(wgcna_dat)

## WGCNA ----
### Prep ----
# Adjust metabolites for basic covariables
wgcna_dat %<>%
  mutate(across(metab_names, 
                ~rstandard(lm(.x~ sex + FASTING + BLCLIN09 + AGEN09)))) %>%
  select(metab_names)

# The following code is based off the WGCNA package tutorial 
# (https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html)
# and https://rstudio-pubs-static.s3.amazonaws.com/90587_39301ea4b9514d20880b4c9cda81ca4f.html
options(stringsAsFactors = FALSE)

# Check for excessive missingness 
# Given metabolite data are imputed, this is just a sense check)
gsg <- goodSamplesGenes(wgcna_dat, verbose = 3)
gsg$allOK #TRUE

### Construct sample network for outlier detection ----
# wgcna_dat <- scale(wgcna_dat[metab_names])
# Transpose data
A <- adjacency(t(wgcna_dat), corFnc = "bicor", type = "distance")
# Network connectivity
k <- as.numeric(apply(A, 2, sum)) - 1
# Standardised connectivity
Z.k <- scale(k) 
summary(Z.k)

# Assign samples at outliying if Z.k < -4 threshold
thresholdZ.k <- -4
sampleTree_outlier <- hclust(as.dist(1 - A), method = "average")

# Plot potential outliers
#pdf("analysis/modules/sample_network.pdf")
#plot(sampleTree_outlier, main = "Sample clustering to detect outliers", sub = "", xlab = "",
#     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
#
#dev.off()

#### Remove outlying samples ----
remove.samples <- Z.k < thresholdZ.k | is.na(Z.k)
table(remove.samples) # 10 = TRUE
wgcna_dat <-  wgcna_dat[!remove.samples,]

# Recompute sample network 
A2 <- adjacency(t(wgcna_dat), corFnc = "bicor", type = "distance")
k2 <- as.numeric(apply(A2, 2, sum)) - 1
Z.k2 <- scale(k2)

datExpr <- wgcna_dat
dim(datExpr)

### Coexpression modules ----
#### Choosing a soft-thresholding power ----
# to meet a scale-free topology  of 0.85, but maximizing mean connectivity
powers <- c(c(1:30))
sft <- pickSoftThreshold(
  datExpr,
  powerVector = powers,
  networkType = "signed",
  corFnc = "bicor",
  corOptions = list(use = 'p', maxPOutliers = 0.1))

sft$powerEstimate # 9

# Plot the results
pdf("analysis/modules/soft-threshold.pdf")
sizeGrWindow(9, 5)
par(mfrow = c(1, 2))
cex1 = 0.9

# Plot scale-free topology fit as a function of the soft-thresholding power
pdf("analysis/modules/soft-threshold.pdf")
plot(ft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Model Fit, signed R^2",
  type = "n",
  main = paste("Scale independence"))

text(sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  labels = powers,
  cex = cex1,
  col = "red")

abline(h = 0.85, col = "red") 

# Plot mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  xlab = "Soft Threshold (power)",
  ylab = "Mean Connectivity",
  type = "n",
  main = paste("Mean connectivity"))

text(sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  labels = powers,
  cex = cex1,
  col = "red")

dev.off()

#### Calculate adjacency matrix using soft thresholding power of 9 ----
softPower <- 9
adjacency <- adjacency(datExpr,
                       power = softPower,
                       type = "signed",
                       corFnc = "bicor")

#### Transform adjacency into topological overlap matrix & calculate dissimilarity ----
TOM <- TOMsimilarity(adjacency, TOMType = "signed")
dissTOM <- 1 - TOM

#### Hierarchical clustering ----
metTree <- hclust(as.dist(dissTOM), method = "average")

# Plot the resulting dendrogram
sizeGrWindow(12, 9)
plot(metTree,
  xlab = "",
  sub = "",
  main = "Metabolite clustering on TOM-based dissimilarity",
  labels = FALSE,
  hang = 0.04)

#### Define modules using dynamic tree cut ----
dynamicMods <- cutreeDynamic(
  dendro = metTree,
  method = "hybrid",
  distM = dissTOM,
  deepSplit = 4,
  minClusterSize = 20)

table(dynamicMods) # 14 modules + grey

# Convert module numbers to (arbitrary) colours
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)

# Plot the dendrogram with colors underneath
pdf("analysis/modules/dendrograms.pdf")

plotDendroAndColors(metTree,
                    dynamicColors,
                    "Dynamic Tree Cut",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05,
                    main = "Metabolite dendrogram and module colors")

dev.off()

#### Calculate module eigenvalues ----
MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs <- MEList$eigengenes

####  Check to see if modules need merging ----
# Calculate module eigenvalue dissimilarity
MEDiss <- 1 - cor(MEs)

# Cluster module eigenvalues
METree <- hclust(as.dist(MEDiss), method = "average")

# Plot
pdf("analysis/modules/merging.pdf")
plot(METree,
     main = "Clustering of module eigenvalues",
     xlab = "",
     sub = "")

MEDissThres <- 0.25 # add cut off line for merging
abline(h = MEDissThres, col = "red") # no modules need merging - all above cut off
dev.off()

### Calculate module membership ----
datKME <- signedKME(datExpr, MEs, outputColumnName = "MM.", corFnc = "bicor")
met_MM <- cbind("metabolon_ID" = row.names(datKME), dynamicColors, datKME) # module membership for each met

# Keep only the module membership for the module the metabolite belongs to
sub_pathway_list <- read.csv("dat/other/sub_pathway_list.csv") %>%
  select(metabolon_ID, BIOCHEMICAL) # for metabolite name

met_MM %<>%
  mutate(MM = lapply(row.names(met_MM),
                      function (x) {
                        col <- met_MM[x, "dynamicColors"]
                        met_MM[x, paste0("MM.", col)]
                      })) %>%
  left_join(sub_pathway_list) %>%
  select(BIOCHEMICAL, module = dynamicColors, MM) %>%
  mutate(MM = as.numeric(MM))

glimpse(met_MM)
write.csv(met_MM, "dat/other/metab_MM.csv", row.names = F)

### Save ----
saveRDS(dynamicColors, "dat/other/dynamicColors.Rdata")
saveRDS(datExpr, "dat/other/datExpr.Rdata")
saveRDS(datKME, "dat/other/datKME.Rdata")
saveRDS(MEs, "dat/other/MEs.Rdata")

## ORA ----
module <- unique(dynamicColors)

# Metabolites & sub/super pathways
sub_pathway_list <- read.csv("dat/other/sub_pathway_list.csv") %>%
  select(metabolon_ID, BIOCHEMICAL, SUB.PATHWAY)

###  Set up ----
# Define pathway reference dataset
ref_met <- sub_pathway_list %>%
  filter(SUB.PATHWAY != "") %>% # remove unknown pathway
  group_by(SUB.PATHWAY) %>%
  filter(n() >= 5) %>% # filter for pathways with 5 or more metabolites
  summarise(n_ref = n()) # count n metabolites per pathway

n_pathways <- nrow(ref_met) # 53 pathways
ref_total_metabolites <- sum(ref_met$n_ref) # 739 metabolites

# Define module pathway data and merge
hyper_dat <- met_MM %>%
  select(BIOCHEMICAL, module) %>%
  left_join(sub_pathway_list, "BIOCHEMICAL") %>%
  group_by(module, SUB.PATHWAY) %>%
  filter(SUB.PATHWAY != "") %>% # remove unknown pathway
  summarise(n_data = n()) %>% # n metabolites per pathway in each module
  inner_join(ref_met, "SUB.PATHWAY")

### Analysis ----
list_of_res <- list()

for (i in module) {
  df <- hyper_dat %>%  # filter data for module
    filter(module == i)
  
  pathway <- as.character(unique(df$SUB.PATHWAY)) 
  
  hyper_results <- matrix(nrow = length(pathway), ncol = 1)
  rownames(hyper_results) <- pathway
  colnames(hyper_results) <- "p_val"
  
  for (j in pathway) {
    
    # number of metabolites in pathway in reference data
    ref_total_pathway <- df$n_ref[which(df$SUB.PATHWAY == j)] 
    # total number of metabolites in module
    sample_total_metabolites <- sum(df$n_data) 
    # number of metabolites in pathway in module
    sample_total_pathway <- df$n_data[which(df$SUB.PATHWAY == j)] 
    
    hyper_results[j, 1] <- as.numeric(phyper(
      sample_total_pathway - 1,
      ref_total_pathway,
      ref_total_metabolites - ref_total_pathway,
      sample_total_metabolites,
      lower.tail = FALSE
    )
    )
  }
  list_of_res[[paste(i)]] <- hyper_results
}

## Save ---
# Supplement
list_of_res <- lapply(list_of_res, data.frame)
write.xlsx(list_of_res, file = "results/ORA_module.xlsx")

# Format for figure
ORA_pathway <- list_of_res %>% 
  do.call(rbind, .) %>%
  mutate("pathway" = row.names(.)) %>% 
  separate(col = "pathway", into = c("module", "pathway"), sep = "\\.") %>%
  mutate(padj = pmin(1, p_val * n_pathways)) %>% 
  select(module, pathway, p_val, padj)

row.names(ORA_pathway) <- NULL
saveRDS(ORA_pathway, "dat/other/ORA_pathway.rds")

# End ----
