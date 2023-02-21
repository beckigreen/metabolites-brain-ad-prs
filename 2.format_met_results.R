#=======================#
#   2.Format results    #
#=======================#

## Set up ----
filepath <- ""
setwd(filepath)

### Load libraries ----
library(dplyr)
library(magrittr)
library(openxlsx)
library(purrr)
library(tidyr)

## Format results (Module) ----
### Brain volume ----
brainvol_mod <- read.csv("results/module/brain_vol_module.csv") %>%
  rename("Module"="X") %>%
  mutate(Module=gsub("ME", "", Module)) 

# FDR correction
brainvol_res_mod <- brainvol_mod %>% 
  select(Module, p1, p2, p3) %>%
  # Stack p-values for FDR correction
  pivot_longer(!Module, names_to="model", values_to="p") %T>% 
  # Check
  {head(.) %>% print} %>%
  # FDR correction
  mutate(padj=p.adjust(p, "fdr")) %>% 
  # Convert back to wide format
  pivot_wider(., names_from="model", values_from=!Module) %>% 
  # Select FDR corrected pvalues
  select(Module, padj_p1, padj_p2, padj_p3) %>% 
  # Include in results df
  merge(brainvol_mod, by="Module") %>% 
  # Reorder
  subset(select=c(1,5:8, 2, 9:12, 3, 13:16, 4)) 

### Hippocampal volume ----
hippvol_mod <- read.csv("results/module/hippo_vol_module.csv") %>%
  rename("Module"="X") %>%
  mutate(Module=gsub("ME", "", Module)) 

# FDR correction
hippvol_res_mod <- hippvol_mod %>% 
  select(Module, p1, p2, p3) %>%
  pivot_longer(!Module, names_to="model", values_to="p") %T>% 
  {head(.) %>% print} %>%
  mutate(padj=p.adjust(p, "fdr")) %>% 
  pivot_wider(., names_from="model", values_from=!Module) %>% 
  select(Module, padj_p1, padj_p2, padj_p3) %>% 
  merge(hippvol_mod, by="Module") %>% 
  subset(select=c(1,5:8, 2, 9:12, 3, 13:16, 4)) 

### Amyloid status  ----
amyloid_mod <- read.csv("results/module/AmyloidStatus_module.csv") %>%
  rename("Module"="X") %>%
  mutate(Module=gsub("ME", "", Module)) 

amyloid_res_mod <- amyloid_mod %>% 
  select(Module, p1, p2, p3) %>%
  pivot_longer(!Module, names_to="model", values_to="p") %T>% 
  {head(.) %>% print} %>%
  mutate(padj=p.adjust(p, "fdr")) %>% 
  pivot_wider(., names_from="model", values_from=!Module) %>% 
  select(Module, padj_p1, padj_p2, padj_p3) %>% 
  merge(amyloid_mod, by="Module") %>% 
  subset(select=c(1,5:8, 2, 9:12, 3, 13:16, 4)) 

### Save ----
list_of_datasets <- list("Brain volume"=brainvol_res_mod, 
                         "Hippocampal volume"=hippvol_res_mod, 
                         "Amyloid Status"=amyloid_res_mod)
write.xlsx(list_of_datasets, file="results/module/insight_modules.xlsx")

## Format results (Hubs) ----
### Read in pathway and module information ----
metab_MM <- read.csv("dat/other/metab_MM.csv") 

# Metabolites & sub/super pathways
pathway <- read.csv("dat/other/sub_pathway_list.csv")

# Hubs from previous paper
cog_hubs <- read.csv("dat/other/hub_id.key.csv") %>%
  pull(BIOCHEMICAL) %>%
  as.character()

### Brain volume ----
brainvol_met <- read.csv("results/metab/brain_vol_metabolites.csv") %>%
  rename("metabolon_ID"="X") %>%
  # Add pathway & module membership info
  left_join(pathway, "metabolon_ID") %>%
  left_join(metab_MM, "BIOCHEMICAL") %>%
  # Indicate whether previous hub
  mutate(`Previous hub`=ifelse(BIOCHEMICAL %in% cog_hubs, 
                                 "Y",
                                 "N")) %>%
  # Reorder
  subset(select=c(16:19, 2:13, 20:21)) 

# FDR correction
brainvol_res_met <- brainvol_met %>% 
  select(BIOCHEMICAL, p1, p2, p3) %>%
  pivot_longer(!BIOCHEMICAL, names_to="model", values_to="p") %T>% 
  {head(.) %>% print} %>%
  mutate(padj=p.adjust(p, "fdr")) %>%
  pivot_wider(., names_from="model", values_from=!BIOCHEMICAL) %>% 
  select(BIOCHEMICAL, padj_p1, padj_p2, padj_p3) %>% 
  merge(brainvol_met, by="BIOCHEMICAL") %>%
  subset(select=c(1,5:11, 2, 12:15, 3, 16:19, 4, 20:21))

### Hippocampal volume ----
hippvol_met <- read.csv("results/metab/hippo_vol_metabolites.csv") %>%
  rename("metabolon_ID"="X") %>%
  left_join(pathway, "metabolon_ID") %>%
  left_join(metab_MM, "BIOCHEMICAL") %>%
  mutate(`Previous hub`=ifelse(BIOCHEMICAL %in% cog_hubs, 
                               "Y",
                               "N")) %>%
  subset(select=c(16:19, 2:13, 20:21)) 

# FDR correction
hippvol_res_met <- hippvol_met %>% 
  select(BIOCHEMICAL, p1, p2, p3) %>%
  pivot_longer(!BIOCHEMICAL, names_to="model", values_to="p") %T>% 
  {head(.) %>% print} %>%
  mutate(padj=p.adjust(p, "fdr")) %>%
  pivot_wider(., names_from="model", values_from=!BIOCHEMICAL) %>% 
  select(BIOCHEMICAL, padj_p1, padj_p2, padj_p3) %>% 
  merge(hippvol_met, by="BIOCHEMICAL") %>%
  subset(select=c(1,5:11, 2, 12:15, 3, 16:19, 4, 20:21))

### Amyloid status ----
amyloid_met <- read.csv("results/metab/AmyloidStatus_metabolites.csv") %>%
  rename("metabolon_ID"="X") %>%
  left_join(pathway, "metabolon_ID") %>%
  left_join(metab_MM, "BIOCHEMICAL") %>%
  mutate(`Previous hub`=ifelse(BIOCHEMICAL %in% cog_hubs, 
                               "Y",
                               "N")) %>%
  subset(select=c(16:19, 2:13, 20:21)) 

# FDR correction
amyloid_res_met <- amyloid_met %>% 
  select(BIOCHEMICAL, p1, p2, p3) %>%
  pivot_longer(!BIOCHEMICAL, names_to="model", values_to="p") %T>% 
  {head(.) %>% print} %>%
  mutate(padj=p.adjust(p, "fdr")) %>%
  pivot_wider(., names_from="model", values_from=!BIOCHEMICAL) %>% 
  select(BIOCHEMICAL, padj_p1, padj_p2, padj_p3) %>% 
  merge(amyloid_met, by="BIOCHEMICAL") %>%
  subset(select=c(1,5:11, 2, 12:15, 3, 16:19, 4, 20:21))

### Save ----
list_of_datasets <- list("Brain volume"=brainvol_res_met, 
                         "Hippocampal volume"=hippvol_res_met, 
                         "Amyloid status"=amyloid_res_met)
write.xlsx(list_of_datasets, file="results/metab/insight_metabolites.xlsx")

### Extract list of hubs ----
pull_hubs <- function(dat) {
  hubs <- dat %>% 
    filter_at(vars(padj_p1, padj_p2, padj_p3), any_vars(.<0.05)) %>%
    pull(BIOCHEMICAL) %>%
    as.character()
  return(hubs)
}

c(pull_hubs(brainvol_res_met),
  pull_hubs(hippvol_res_met),
  pull_hubs(amyloid_res_met)) %>%
  unique() %T>%
  {length(.) %>% print} %>%
  write.table("dat/other/FDRhub.csv", row.names = F)
  
# End ----