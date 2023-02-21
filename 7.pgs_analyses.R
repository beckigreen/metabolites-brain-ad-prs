#====================#
#   7.PGS Analyses   #
#====================#

## Set up ----
filepath <- ""
setwd(filepath)

### Load libraries ----
library(dplyr)
library(magrittr)
library(purrr)
library(tidyr)

### Load data ----
# Metabolite & outcome imputed data
dat <- readRDS("dat/ins_met_nodementia.Rdata")

# Metabolites & sub/super pathways
pathway <- read.csv("dat/other/sub_pathway_list.csv")

# PGS results 
pgs <- read.table("analysis/pgs/NSHD_PRS.all_score", header =T) %>% 
  rename("nshdid" = "IID", "PGS_APOE" = "Pt_5e.08", "PGS_APOE0.1" = "Pt_0.1")
pgsnoapoe <- read.table("analysis/pgs/NSHD_PRS_noAPOE.all_score", header =T) %>% 
  rename("nshdid" = "IID", "PGS_noAPOE" = "Pt_5e.08", "PGS_noAPOE0.1" = "Pt_0.1")

PCs <- read.csv("dat/genetics/target/pcs.csv")

### Merge data ----
pgs_dat <- reduce(list(dat$data, pgs, pgsnoapoe, PCs), left_join, "nshdid")

pgs_dat %>% 
  filter(!(is.na(PGS_APOE) & is.na(PGS_noAPOE))) %T>%
  {nrow(.) %>% print}

### Analysis (module) ----
module <- c("MEyellow", "MEblue", "MEbrown")

module_results <- matrix(data = NA, nrow = 3, ncol =16)
rownames(module_results) <- module
colnames(module_results) <- c("beta_APOE", "L95_APOE", "U95_APOE", "p_APOE",
                              "beta_APOE0.1", "L95_APOE0.1", "U95_APOE0.1", "p_APOE0.1",
                              "beta_noAPOE", "L95_noAPOE", "U95_noAPOE", "p_noAPOE",
                              "beta_noAPOE0.1", "L95_noAPOE0.1", "U95_noAPOE0.1", "p_noAPOE0.1")
for (i in module){
  print(i)

  print("APOE (pT=5x10-8)...")
  frmla_apoe <- paste0("scale(", i, ")", 
  " ~ scale(PGS_APOE) + scale(PC1) + scale(PC2) + scale(PC3) + 
  scale(PC4) + scale(PC5) + scale(PC6) + scale(PC7)") 
  modelFit_apoe <- lm(as.formula(frmla_apoe), data = pgs_dat)
  
  res_apoe <- summary(modelFit_apoe, conf.int = TRUE)
  beta_apoe <- res_apoe$coefficients[2,1]
  L95_apoe <- confint(modelFit_apoe)[2,1]
  U95_apoe <- confint(modelFit_apoe)[2,2]
  p_apoe <- res_apoe$coefficients[2,4]
  
  module_results[i,1] <- beta_apoe
  module_results[i,2] <- L95_apoe
  module_results[i,3] <- U95_apoe
  module_results[i,4] <- p_apoe
  
  print("APOE (pT=0.1)...")
  frmla_apoe0.1 <- paste0("scale(", i, ")", 
  " ~ scale(PGS_APOE0.1) + scale(PC1) + scale(PC2) + scale(PC3) + 
  scale(PC4) + scale(PC5) + scale(PC6) + scale(PC7)")
  modelFit_apoe0.1 <- lm(as.formula(frmla_apoe0.1), data = pgs_dat)
  
  res_apoe0.1 <- summary(modelFit_apoe0.1, conf.int = TRUE)
  beta_apoe0.1 <- res_apoe0.1$coefficients[2,1]
  L95_apoe0.1 <- confint(modelFit_apoe0.1)[2,1]
  U95_apoe0.1 <- confint(modelFit_apoe0.1)[2,2]
  p_apoe0.1 <- res_apoe0.1$coefficients[2,4]
  
  module_results[i,5] <- beta_apoe0.1
  module_results[i,6] <- L95_apoe0.1
  module_results[i,7] <- U95_apoe0.1
  module_results[i,8] <- p_apoe0.1
  
  print("no APOE (pT=5x10-8)...")
  frmla_noapoe <- paste0("scale(", i, ")", 
  " ~ scale(PGS_noAPOE) + scale(PC1) + scale(PC2) + scale(PC3) + 
  scale(PC4) + scale(PC5) + scale(PC6) + scale(PC7)") 
  modelFit_noapoe <- lm(as.formula(frmla_noapoe), data = pgs_dat)
  
  res_noapoe <- summary(modelFit_noapoe, conf.int = TRUE)
  beta_noapoe <- res_noapoe$coefficients[2,1]
  L95_noapoe <- confint(modelFit_noapoe)[2,1]
  U95_noapoe <- confint(modelFit_noapoe)[2,2]
  p_noapoe <- res_noapoe$coefficients[2,4]
  
  module_results[i,9] <- beta_noapoe
  module_results[i,10] <- L95_noapoe
  module_results[i,11] <- U95_noapoe
  module_results[i,12] <- p_noapoe
  
  print("no APOE (pT=0.1)...")
  frmla_noapoe0.1 <- paste0("scale(", i, ")", 
  " ~ scale(PGS_noAPOE0.1) + scale(PC1) + scale(PC2) + scale(PC3) + 
  scale(PC4) + scale(PC5) + scale(PC6) + scale(PC7)") 
  modelFit_noapoe0.1 <- lm(as.formula(frmla_noapoe0.1), data = pgs_dat)
  
  res_noapoe0.1 <- summary(modelFit_noapoe0.1, conf.int = TRUE)
  beta_noapoe0.1 <- res_noapoe0.1$coefficients[2,1]
  L95_noapoe0.1 <- confint(modelFit_noapoe0.1)[2,1]
  U95_noapoe0.1 <- confint(modelFit_noapoe0.1)[2,2]
  p_noapoe0.1 <- res_noapoe0.1$coefficients[2,4]
  
  module_results[i,13] <- beta_noapoe0.1
  module_results[i,14] <- L95_noapoe0.1
  module_results[i,15] <- U95_noapoe0.1
  module_results[i,16] <- p_noapoe0.1
  
}

### FDR correction ----
module_results_adj <- as.data.frame(module_results) %>% 
  mutate(Module = row.names(.))  %>% 
  pivot_longer(c("p_APOE", "p_noAPOE", "p_APOE0.1", "p_noAPOE0.1"), 
               names_to = "analysis", values_to = "p") %>% # stack p values
  mutate(padj = p.adjust(p, "fdr")) %>%  # fdr correction
  pivot_wider(names_from = "analysis", values_from = c(p, padj)) %>% 
  select(c(13, 1:3, 14, 18, 4:6, 16, 20, 7:9, 15, 19, 10:12, 17, 21)) #reorder

names(module_results_adj) <- gsub(names(module_results_adj),pattern= "_p", replacement="")

### Save ----
write.csv(module_results_adj, "results/pgs/PGS_modules.csv", row.names=F)

### Analysis (hubs) ----
hubs <- read.csv("dat/other/FDRhub.csv") %>%
  rename(BIOCHEMICAL = x) %>%
  as.data.frame() %>%
  left_join(., pathway) %>%
  select(BIOCHEMICAL, metabolon_ID)

metab_results <- matrix(data = NA, nrow = length(hubs$metabolon_ID), ncol =16)
rownames(metab_results) <- hubs$metabolon_ID
colnames(metab_results) <- c("beta_APOE", "L95_APOE", "U95_APOE", "p_APOE",
                             "beta_APOE0.1", "L95_APOE0.1", "U95_APOE0.1", "p_APOE0.1",
                             "beta_noAPOE", "L95_noAPOE", "U95_noAPOE", "p_noAPOE",
                             "beta_noAPOE0.1", "L95_noAPOE0.1", "U95_noAPOE0.1", "p_noAPOE0.1")
for (i in hubs$metabolon_ID){
  metab_name <- hubs$BIOCHEMICAL[which(hubs$metabolon_ID == i)]
  print(metab_name)
  
  print("APOE (pT=5x10-8)...")
  frmla_apoe <- paste0("scale(", i, ")", 
                       " ~ scale(PGS_APOE) + scale(AGEN09) + sex + FASTING + 
                       BLCLIN09 + scale(PC1) + scale(PC2) + scale(PC3) + 
                       scale(PC4) + scale(PC5) + scale(PC6) + scale(PC7)") 
  modelFit_apoe <- lm(as.formula(frmla_apoe), data = pgs_dat)
  
  res_apoe <- summary(modelFit_apoe, conf.int = TRUE)
  beta_apoe <- res_apoe$coefficients[2,1]
  L95_apoe <- confint(modelFit_apoe)[2,1]
  U95_apoe <- confint(modelFit_apoe)[2,2]
  p_apoe <- res_apoe$coefficients[2,4]
  
  metab_results[i,1] <- beta_apoe
  metab_results[i,2] <- L95_apoe
  metab_results[i,3] <- U95_apoe
  metab_results[i,4] <- p_apoe
  
  print("APOE (pT=0.1)...")
  frmla_apoe0.1 <- paste0("scale(", i, ")", 
                          " ~ scale(PGS_APOE0.1) + scale(AGEN09) + sex + FASTING + 
                       BLCLIN09 + scale(PC1) + scale(PC2) + scale(PC3) + 
                       scale(PC4) + scale(PC5) + scale(PC6) + scale(PC7)")
  modelFit_apoe0.1 <- lm(as.formula(frmla_apoe0.1), data = pgs_dat)
  
  res_apoe0.1 <- summary(modelFit_apoe0.1, conf.int = TRUE)
  beta_apoe0.1 <- res_apoe0.1$coefficients[2,1]
  L95_apoe0.1 <- confint(modelFit_apoe0.1)[2,1]
  U95_apoe0.1 <- confint(modelFit_apoe0.1)[2,2]
  p_apoe0.1 <- res_apoe0.1$coefficients[2,4]
  
  metab_results[i,5] <- beta_apoe0.1
  metab_results[i,6] <- L95_apoe0.1
  metab_results[i,7] <- U95_apoe0.1
  metab_results[i,8] <- p_apoe0.1
  
  print("no APOE (pT=5x10-8)...")
  frmla_noapoe <- paste0("scale(", i, ")", 
                         " ~ scale(PGS_noAPOE) + scale(AGEN09) + sex + FASTING + 
                         BLCLIN09 + scale(PC1) + scale(PC2) + scale(PC3) + 
                         scale(PC4) + scale(PC5) + scale(PC6) + scale(PC7)")
  modelFit_noapoe <- lm(as.formula(frmla_noapoe), data = pgs_dat)
  
  res_noapoe <- summary(modelFit_noapoe, conf.int = TRUE)
  beta_noapoe <- res_noapoe$coefficients[2,1]
  L95_noapoe <- confint(modelFit_noapoe)[2,1]
  U95_noapoe <- confint(modelFit_noapoe)[2,2]
  p_noapoe <- res_noapoe$coefficients[2,4]
  
  metab_results[i,9] <- beta_noapoe
  metab_results[i,10] <- L95_noapoe
  metab_results[i,11] <- U95_noapoe
  metab_results[i,12] <- p_noapoe
  
  print("no APOE (pT=0.1)...")
  frmla_noapoe0.1 <- paste0("scale(", i, ")", 
                            " ~ scale(PGS_noAPOE0.1) + scale(AGEN09) + sex + FASTING + 
                            BLCLIN09 + scale(PC1) + scale(PC2) + scale(PC3) + 
                            scale(PC4) + scale(PC5) + scale(PC6) + scale(PC7)") 
  modelFit_noapoe0.1 <- lm(as.formula(frmla_noapoe0.1), data = pgs_dat)
  
  res_noapoe0.1 <- summary(modelFit_noapoe0.1, conf.int = TRUE)
  beta_noapoe0.1 <- res_noapoe0.1$coefficients[2,1]
  L95_noapoe0.1 <- confint(modelFit_noapoe0.1)[2,1]
  U95_noapoe0.1 <- confint(modelFit_noapoe0.1)[2,2]
  p_noapoe0.1 <- res_noapoe0.1$coefficients[2,4]
  
  metab_results[i,13] <- beta_noapoe0.1
  metab_results[i,14] <- L95_noapoe0.1
  metab_results[i,15] <- U95_noapoe0.1
  metab_results[i,16] <- p_noapoe0.1
  
  # Rownames to metabolite name rather than ID
  rownames(metab_results)[rownames(metab_results) == i] <- as.character(metab_name)
  
}

### FDR correction ----
metab_results_adj <- as.data.frame(metab_results) %>% 
  mutate(BIOCHEMICAL = row.names(.))  %>% 
  pivot_longer(c("p_APOE", "p_noAPOE", "p_APOE0.1", "p_noAPOE0.1"), 
               names_to = "analysis", values_to = "p") %>% # stack p values
  mutate(padj = p.adjust(p, "fdr")) %>%  # fdr correction
  pivot_wider(names_from = "analysis", values_from = c(p, padj)) %>% 
  select(c(13, 1:3, 14, 18, 4:6, 16, 20, 7:9, 15, 19, 10:12, 17, 21)) #reorder

names(metab_results_adj) <- gsub(names(metab_results_adj),pattern= "_p", replacement="")

### Save ----
write.csv(metab_results_adj, "results/pgs/PGS_metabs.csv", row.names=F)

# End ----