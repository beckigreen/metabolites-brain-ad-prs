#==========================#
#  3. Base (& target) QC   #
#==========================#

## Set up ----
filepath <- ""
setwd(filepath)

### Load libraries ----
library(data.table)
library(dplyr)
library(magrittr)
library(purrr)

### Load data ----
base <- fread("dat/genetics/base/Kunkle_etal_Stage1_results.txt", header = T)
colnames(base) <- c("CHR", "BP", "SNP", "A1", "A2", "BETA", "SE", "P")
glimpse(base)

## Checks (base data) ----
# NAs & duplicates
map(base, ~sum(is.na(.))) # some missing in SNP col
sum(duplicated(base$SNP)) # some duplicates in SNP col

base %>% 
  filter(duplicated(SNP)) %>% 
  select(SNP) %>%
  table(useNA = "always") # all duplicates are NAs

# RSIDs
base %>%
  filter(!grepl("rs", SNP)) %>%
  nrow() # some SNPs with no rsid

# Alleles
base %>% 
  filter(nchar(A1) > 1 | nchar(A2) > 1) %>%
  nrow() # some multichar alleles 

base %>%
  filter(nchar(A1) == 1 & nchar(A2) == 1) %>%
  select(A1, A2) %>%
  map(table) # all looks fine

### QC base data ----
base_qc <- base %>% 
  filter(!is.na(SNP)) %>% 
  filter(nchar(A1) == 1 & nchar(A2) == 1) %>%
  filter(grepl("rs", SNP)) %>%
  mutate_at(c("CHR", "BP", "A2", "A2"), as.character) %>%
  mutate_at(c("BETA", "SE", "P"), as.numeric)

### Problem --- P vals are converted to zero as very small
p0snps <- base_qc %>%
  filter(P == 0) %T>%
  {nrow(.) %>% print} %>% # 18 SNPs
  pull(SNP)
  
base %>% 
  filter(SNP %in% p0snps) %>%
  mutate(Z=BETA/SE) %>% 
  arrange(-abs(Z)) %>%
  select(SNP, Z, P) 

# very small P values, smallest for rs429358 (1.17e-881)
# highest absolute Z score (63.6)

write.table(p0snps, "dat/genetics/base/p0snps.txt", col.names=F, row.names=F, quote=F)

# This problem also lies if left as character and read in by PRSice
# Issue is that PRSice will then be unstable with..
# ...which of the SNPs selected for PGS (if all are in LD)
# Thanks to Dr Sam Choi for their help solving this!

# Run the following using plink 
# plink \
#--bfile dat/genetics/target/NSHD_QCed \
#--extract dat/genetics/base/p0snps.txt \
#--r2 \
#--out dat/genetics/target/ldmat

# All p=0 SNPs in LD
# SNP with highest absolute Z-score / lowest p-value is rs429358
# others can be converted to NA
base_qc$P[which(base_qc$P == 0)] <- NA
base_qc$P[which(base_qc$SNP == "rs429358")] <- 0 # NB: SNP is in target

## Save ----
write.table(base_qc, 
            file = "dat/genetics/base/Kunkle_Stage1_post-qc_semi_clumped.txt", 
            quote = F, row.names = F, col.names = T)

## Checks (target data) ----
target <- fread("dat/genetics/target/NSHD_QCed.bim")
colnames(target) <- c("CHR", "SNP", "CM", "BP", "A1", "A2")

map(target, ~sum(is.na(.)))
sum(duplicated(target$SNP))

target %>%
  filter(!grepl("rs", SNP)) %>%
  nrow() 

target %>% 
  filter(nchar(A1) > 1 | nchar(A2) > 1) %>%
  nrow() 

target %>%
  select(A1, A2) %>%
  map(table) 

map(target, class)

# All looks fine

# Pull SNPs in APOE region
target %>% 
  filter(CHR == 19) %>% 
  filter(BP >= 44912079 & BP <= 45912079) %>% 
  pull(SNP) %T>%
  {head(.) %>% print} %>%
  write.table("dat/genetics/target/APOE_snps_NSHD.txt", col.names=F, row.names=F, quote=F) 

## End ----