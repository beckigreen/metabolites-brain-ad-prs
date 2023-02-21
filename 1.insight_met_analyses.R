#============================#
#  1.Analyses: Module & Hub  #
#============================#

## Set up ----
filepath <- ""
setwd(filepath)

### Load libraries ----
library(dplyr)
library(magrittr)
library(mice)

### Read in data ----
# Metabolite & outcome imputed data
dat <- readRDS("dat/ins_met_nodementia.Rdata")

# Metabolite module memberships
metab_MM <- read.csv("dat/other/metab_MM.csv")

# Metabolites & sub/super pathways
pathway <- read.csv("dat/other/sub_pathway_list.csv")

# Hubs from previous paper
cog_hubs <- read.csv("dat/other/hub_id.key.csv") %>%
  pull(metabolon_ID) %>%
  as.character()

### Functions ----
# Regression diagnostics with mids object
mids_diagnostics <- function(modelFit, predictor_name="") {
  res <- sapply(modelFit$analyses, residuals)
  res <- apply(res, 1, mean)
  fit <- sapply(modelFit$analyses, fitted)
  fit <- apply(fit, 1, mean)
  res_standard <- sapply(modelFit$analyses, rstandard)
  res_standard <- apply(res_standard, 1, mean)
  
  plot(fit,
       res,
       main=paste(predictor_name),
       font.sub=2)
  qqnorm(res_standard,
         main=paste(predictor_name),
         font.sub=2)
  qqline(res_standard)
}

## Module analysis -----
### Brain volumes ----
vol_outcome <- c("brain_vol", "hippo_vol")
module <-c(
    "MEyellow",
    "MEcyan",
    "MEbrown",
    "MEturquoise",
    "MEpurple",
    "MEmagenta",
    "MEred",
    "MEblue",
    "MEgreenyellow",
    "MEsalmon",
    "MEpink",
    "MEtan",
    "MEblack",
    "MEgreen"
  )

# Nested models
# NB: modules already adjusted for fasting, age at blood collection & clinic
# Model 1=APOE4, intracranial volume, age at scan
mod1_volume <- " + scale(spm_tiv_vol) + apoe4 + scale(age_scan_vol)"

# Model 2=model 1 + BMI + lipid medication
mod2_volume <- paste0(mod1_volume, "+ scale(bmi09) + MEDChol_63x")

# Model 3=model 2 + childhood cognition + educational attainment + child and adult SEP
mod3_volume <- paste0(mod2_volume, "+ scale(cog15h) + LHQR + chsc + SC1553")

# Analysis
for (i in vol_outcome) {
  print(i)
  pdf(paste0("results/module/", i,"_module_lm.pdf")) # set up file for regression diagnostics plots
  module_results <- matrix(data=NA, nrow=14, ncol=12) # create results matrix
  rownames(module_results) <- module
  colnames(module_results) <- c(
    # model 1
    "beta1",
    "L951",
    "U951",
    "p1",
    # model 2
    "beta2",
    "L952",
    "U952",
    "p2",
    # model 3
    "beta3",
    "L953",
    "U953",
    "p3"
  )
  
  for (j in module) {
    print(j)
    
    print("Model 1...")
    frmla <-paste0("scale(", i,")", " ~ scale(", j,")", mod1_volume)
    modelFit_m1 <- with(dat, lm(as.formula(frmla)))
    
    res_m1 <- summary(pool(modelFit_m1), conf.int=TRUE)
    beta_m1 <- res_m1[2, 2]
    L95_m1 <- res_m1[2, 7]
    U95_m1 <- res_m1[2, 8]
    p_m1 <- res_m1[2, 6]
    
    module_results[j,1] <- beta_m1
    module_results[j,2] <- L95_m1
    module_results[j,3] <- U95_m1
    module_results[j,4] <- p_m1
    
    mids_diagnostics(modelFit_m1, predictor_name=paste(j, "- Model 1")) # regression diagnostics
    
    print("Model 2...")
    frmla2 <-paste0("scale(", i,")", " ~ scale(", j,")", mod2_volume)
    modelFit_m2 <- with(dat, lm(as.formula(frmla2)))
    
    res_m2 <- summary(pool(modelFit_m2), conf.int=TRUE)
    beta_m2 <- res_m2[2, 2]
    L95_m2 <- res_m2[2, 7]
    U95_m2 <- res_m2[2, 8]
    p_m2 <- res_m2[2, 6]
    
    module_results[j,5] <- beta_m2
    module_results[j,6] <- L95_m2
    module_results[j,7] <- U95_m2
    module_results[j,8] <- p_m2
    
    mids_diagnostics(modelFit_m2, predictor_name=paste(j, "- Model 2"))
    
    print("Model 3...")
    frmla3 <- paste0("scale(", i,")", " ~ scale(", j, ")", mod3_volume)
    modelFit_m3 <- with(dat, lm(as.formula(frmla3)))
    
    res_m3 <- summary(pool(modelFit_m3), conf.int=TRUE)
    beta_m3 <- res_m3[2, 2]
    L95_m3 <- res_m3[2, 7]
    U95_m3 <- res_m3[2, 8]
    p_m3 <- res_m3[2, 6]
    
    module_results[j,9] <- beta_m3
    module_results[j,10] <- L95_m3
    module_results[j,11] <- U95_m3
    module_results[j,12] <- p_m3
    
    mids_diagnostics(modelFit_m3, predictor_name=paste(j,"- Model 3"))
    
  }
  print("Saving...")
  write.csv(module_results,
    paste0("results/module/", i,"_module.csv"))
  dev.off()
}

### Amyloid status ----
mod1_amyloid <- " + apoe4 + scale(age_scan_amy)"
mod2_amyloid <- paste0(mod1_amyloid," + scale(bmi09) + MEDChol_63x")
mod3_amyloid <- paste0(mod2_amyloid, " + scale(cog15h) + LHQR + chsc + SC1553")

module_results <- matrix(data=NA, nrow=14, ncol=12)
rownames(module_results) <- module
colnames(module_results) <- c(
  "OR1",
  "L951",
  "U951",
  "p1",
  "OR2",
  "L952",
  "U952",
  "p2",
  "OR3",
  "L953",
  "U953",
  "p3"
)

for (i in module) {
  print(i)
  
  print("Model 1...")
  frmla <- paste0("AmyloidStatus ~ scale(", i,")", mod1_amyloid)
  modelFit_m1 <- with(dat, glm(as.formula(frmla), family="binomial"))
  
  res_m1 <- summary(pool(modelFit_m1), conf.int=TRUE, exponentiate=TRUE)
  OR_m1 <- res_m1[2, 2]
  L95_m1 <- res_m1[2, 7]
  U95_m1 <- res_m1[2, 8]
  p_m1 <- res_m1[2, 6]
  
  module_results[i,1] <- OR_m1
  module_results[i,2] <- L95_m1
  module_results[i,3] <- U95_m1
  module_results[i,4] <- p_m1
  
  print("Model 2...")
  frmla2 <- paste0("AmyloidStatus ~ scale(", i,")", mod2_amyloid)
  modelFit_m2 <- with(dat, glm(as.formula(frmla2), family="binomial"))
  
  res_m2 <- summary(pool(modelFit_m2), conf.int=TRUE, exponentiate=TRUE)
  OR_m2 <- res_m2[2, 2]
  L95_m2 <- res_m2[2, 7]
  U95_m2 <- res_m2[2, 8]
  p_m2 <- res_m2[2, 6]
  
  module_results[i,5] <- OR_m2
  module_results[i,6] <- L95_m2
  module_results[i,7] <- U95_m2
  module_results[i,8] <- p_m2
  
  print("Model 3...")
  frmla3 <- paste0("AmyloidStatus ~ scale(", i,")", mod3_amyloid)
  modelFit_m3 <- with(dat, glm(as.formula(frmla3), family="binomial"))
  
  res_m3 <- summary(pool(modelFit_m3), conf.int=TRUE, exponentiate=TRUE)
  OR_m3 <- res_m3[2, 2]
  L95_m3 <- res_m3[2, 7]
  U95_m3 <- res_m3[2, 8]
  p_m3 <- res_m3[2, 6]
  
  module_results[i,9] <- OR_m3
  module_results[i,10] <- L95_m3
  module_results[i,11] <- U95_m3
  module_results[i,12] <- p_m3
  
}
print("Saving...")
write.csv(module_results, "results/module/AmyloidStatus_module.csv")

## Metabolite analysis ----
### Define hubs ----
hubs <- metab_MM  %>%
  # Filter metabolites to those in modules showing assoc with an outcome
  filter(module %in% c("yellow", "brown", "blue")) %>%
  # Filter for module membership > 0.65
  filter(MM > 0.65) %>%
  # Get ID
  select(BIOCHEMICAL) %>%
  left_join(pathway) %>%
  pull(metabolon_ID) %>%
  unique() %>%
  as.character() %T>%
  {length(.) %>% print}

# Combine with cognition hubs
hubs %<>%
  c(hubs, cog_hubs) %>%
  unique() %T>%
  {length(.) %>% print}

### Brain volumes ----
# Nested models
# Model 1=blood clinic info, APOE4, intracranial volume, age at scan
met1_volume <- "+ scale(AGEN09) + FASTING + BLCLIN09 + sex  + scale(spm_tiv_vol) + apoe4 + scale(age_scan_vol)"

# Model 2=model 1 + BMI + lipid medication
met2_volume  <- paste0(met1_volume, "+ scale(bmi09) + MEDChol_63x")

# Model 3=model 2 + childhood cognition + educational attainment + child and adult SEP
met3_volume  <- paste0(met2_volume, "+ scale(cog15h) + LHQR + chsc + SC1553")

vol_outcome <- c("hippo_vol", "brain_vol")

for (i in vol_outcome) {
  print(i)
  pdf(paste0("results/metab/", i,"_metabolite_lm.pdf"))
  metab_results <- matrix(data=NA,nrow=length(hubs), ncol=12) # create results matrix
  rownames(metab_results) <- hubs
  colnames(metab_results) <- c(
    # model 1
    "beta1",
    "L951",
    "U951",
    "p1",
    # model 2
    "beta2",
    "L952",
    "U952",
    "p2",
    # model 3
    "beta3",
    "L953",
    "U953",
    "p3"
  )
  
  for (j in hubs) {
    print(j)
    
    print("Model 1...")
    frmla <- paste0("scale(", i,")", " ~ scale(", j,")", met1_volume)
    modelFit_m1 <- with(dat, lm(as.formula(frmla)))
    
    res_m1 <- summary(pool(modelFit_m1), conf.int=TRUE)
    beta_m1 <- res_m1[2, 2]
    L95_m1 <- res_m1[2, 7]
    U95_m1 <- res_m1[2, 8]
    p_m1 <- res_m1[2, 6]
    
    metab_results[j,1] <- beta_m1
    metab_results[j,2] <- L95_m1
    metab_results[j,3] <- U95_m1
    metab_results[j,4] <- p_m1
    
    mids_diagnostics(modelFit_m1, predictor_name=paste(j,"- Model 1")) # regression diagnostics
    
    print("Model 2...")
    frmla2 <- paste0("scale(", i,")", " ~ scale(", j,")", met2_volume)
    modelFit_m2 <- with(dat, lm(as.formula(frmla2)))
    
    res_m2 <- summary(pool(modelFit_m2), conf.int=TRUE)
    beta_m2 <- res_m2[2, 2]
    L95_m2 <- res_m2[2, 7]
    U95_m2 <- res_m2[2, 8]
    p_m2 <- res_m2[2, 6]
    
    metab_results[j,5] <- beta_m2
    metab_results[j,6] <- L95_m2
    metab_results[j,7] <- U95_m2
    metab_results[j,8] <- p_m2
    
    mids_diagnostics(modelFit_m2, predictor_name=paste(j,"- Model 2"))
    
    print("Model 3...")
    frmla3 <- paste0("scale(", i,")", " ~ scale(", j,")", met3_volume)
    modelFit_m3 <- with(dat, lm(as.formula(frmla3)))
    
    res_m3 <- summary(pool(modelFit_m3), conf.int=TRUE)
    beta_m3 <- res_m3[2, 2]
    L95_m3 <- res_m3[2, 7]
    U95_m3 <- res_m3[2, 8]
    p_m3 <- res_m3[2, 6]
    
    metab_results[j,9] <- beta_m3
    metab_results[j,10] <- L95_m3
    metab_results[j,11] <- U95_m3
    metab_results[j,12] <- p_m3
    
    mids_diagnostics(modelFit_m3, predictor_name=paste(j,"- Model 3"))
    
  }
  print("Saving...")
  write.csv(metab_results, paste0("results/metab/", i,"_metabolites.csv"))
  dev.off()
}

### Amyloid status ----
met1_amyloid <- "+ scale(AGEN09) + FASTING + BLCLIN09 + sex + apoe4 + scale(age_scan_amy)"
met2_amyloid  <- paste0(met1_amyloid, " + scale(bmi09) + MEDChol_63x")
met3_amyloid  <- paste0(met2_amyloid, " + scale(cog15h) + LHQR + chsc + SC1553")

metab_results <- matrix(data=NA, nrow=length(hubs), ncol=12)
rownames(metab_results) <- hubs
colnames(metab_results) <- c(
  "OR1",
  "L951",
  "U951",
  "p1",
  "OR2",
  "L952",
  "U952",
  "p2",
  "OR3",
  "L953",
  "U953",
  "p3"
)

for (i in hubs) {
  print(i)
  
  print("Model 1...")
  frmla <- paste0("AmyloidStatus ~ scale(", i,")", met1_amyloid)
  modelFit_m1 <- with(dat, glm(as.formula(frmla), family="binomial"))
  
  res_m1 <- summary(pool(modelFit_m1), conf.int=TRUE, exponentiate=TRUE)
  OR_m1 <- res_m1[2, 2]
  L95_m1 <- res_m1[2, 7]
  U95_m1 <- res_m1[2, 8]
  p_m1 <- res_m1[2, 6]
  
  metab_results[i,1] <- OR_m1
  metab_results[i,2] <- L95_m1
  metab_results[i,3] <- U95_m1
  metab_results[i,4] <- p_m1
  
  print("Model 2...")
  frmla2 <- paste0("AmyloidStatus ~ scale(", i,")", met2_amyloid)
  modelFit_m2 <- with(dat, glm(as.formula(frmla2), family="binomial"))
  
  res_m2 <- summary(pool(modelFit_m2), conf.int=TRUE, exponentiate=TRUE)
  OR_m2 <- res_m2[2, 2]
  L95_m2 <- res_m2[2, 7]
  U95_m2 <- res_m2[2, 8]
  p_m2 <- res_m2[2, 6]
  
  metab_results[i,5] <- OR_m2
  metab_results[i,6] <- L95_m2
  metab_results[i,7] <- U95_m2
  metab_results[i,8] <- p_m2
  
  print("Model 3...")
  frmla3 <- paste0("AmyloidStatus ~ scale(", i,")", met3_amyloid)
  modelFit_m3 <- with(dat, glm(as.formula(frmla3), family="binomial"))
 
  res_m3 <- summary(pool(modelFit_m3), conf.int=TRUE, exponentiate=TRUE)
  OR_m3 <- res_m3[2, 2]
  L95_m3 <- res_m3[2, 7]
  U95_m3 <- res_m3[2, 8]
  p_m3 <- res_m3[2, 6]
  
  metab_results[i,9] <- OR_m3
  metab_results[i,10] <- L95_m3
  metab_results[i,11] <- U95_m3
  metab_results[i,12] <- p_m3
  
}
print("Saving...")
write.csv(metab_results, "results/metab/AmyloidStatus_metabolites.csv")

# End ----