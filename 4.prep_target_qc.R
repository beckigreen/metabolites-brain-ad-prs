#=========================#
# 4.Prep target QC files  #
#=========================#

## Set up ----
filepath <- ""
setwd(filepath)

### Load libraries ----
library(dplyr)
library(magrittr)

## Sex checks (checking whether genetic and self-report match) ----
# Genetic QC files performed by Marzia Scelsi
# Two different genotyping chips; "PROBLEM" indicates potential sample mixup
sexcheck_neuro <- read.csv("dat/genetics/target/scrambled_1946_NeuroX2_batch1and2_noCtrls_sex.sexcheck.csv")
sexcheck_drugdev <- read.csv("dat/genetics/target/scrambled_1946_DrugDev_all_noCtrls.sexcheck.csv")

rbind(sexcheck_neuro, sexcheck_drugdev) %>%
  select(FID, IID, PEDSEX) %>%
  distinct() %>%
  write.table("dat/genetics/target/ID_sex.txt",
              sep="\t", 
              row.names = F, 
              col.names = F)

sexcheck_neuro %<>% 
  filter(STATUS == "PROBLEM") %>% 
  select(FID, IID) %T>%
  {nrow(.) %>% print}

sexcheck_drugdev %<>% 
  filter(STATUS == "PROBLEM") %>% 
  select(FID, IID) %T>%
  {nrow(.) %>% print}

sexcheck_all <- rbind(sexcheck_neuro, sexcheck_drugdev) %>% 
  distinct() %T>%
  {nrow(.) %>% print}

# Save
write.table(sexcheck_neuro, 
            "dat/genetics/target/sexcheck_neuro.txt",
            sep="\t", 
            row.names = F, 
            col.names = F)
write.table(sexcheck_drugdev, 
            "dat/genetics/target/sexcheck_drugdev.txt",
            sep="\t", 
            row.names = F, 
            col.names = F)

## Relatedness ----
# Genetic QC files performed by Marzia Scelsi
# CSV lists unrelated IDs (PIHAT <0.1)
rel <- read.csv("dat/genetics/target/scrambled_1946_relatedness0.1.rel.id.csv")
colnames(rel) <- c("FID", "IID")

## Relatedness checks
# See N will be removing overall
rel %T>%
  {nrow(.) %>% print} %>%
  filter(!(FID %in% sexcheck_all$FID)) %>%
  n_distinct() 

# Save as txt
write.table(rel, 
            "dat/genetics/target/rel.txt",
            sep="\t", 
            row.names = F, 
            col.names = F)

## PCs
pcs <- read.csv("dat/genetics/target/scrambled_1946_NeuroX2_DrugDev_all_autsms_AFdiff10-updated-REFfixed_QC_CEU.eigenvec.csv")
colnames(pcs) <- c("nshdid", "IID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")

glimpse(pcs)
write.csv(pcs, "dat/genetics/target/pcs.csv", row.names=F)

# End ----
