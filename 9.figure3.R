#====================#
#      Figure 3      #
#====================#

## Set up ----
filepath <- ""
setwd(filepath)

### Load libraries ----
library(cowplot)
library(dplyr)
library(ggh4x)
library(ggplot2)
library(grid)
library(ggtext)
library(magrittr)
library(openxlsx)
library(purrr)
library(tidyr)

### Read in files ----
# Metabolite module memberships
metab_MM <- read.csv("dat/other/metab_MM.csv")

# Metabolites & sub/super pathways
pathway <- read.csv("dat/other/sub_pathway_list.csv")

# Hubs from previous paper
cog_hubs <- read.csv("dat/other/hub_id.key.csv") %>%
  pull(BIOCHEMICAL) %>%
  as.character()

## Format hub results for plotting ----
bv <- read.xlsx("results/metab/insight_metabolites.xlsx", sheet=1) %>%
  # stack result
  gather(variable, value, 
         c("beta1", "beta2", "beta3", 
           "L951", "L952", "L953",
           "U951", "U952", "U953",
           "p1", "p2", "p3",
           "padj_p1", "padj_p2", "padj_p3")) %>%
  # create col indicating model
  mutate(model = case_when(
      variable == "beta1" | variable ==  "U951"
      | variable == "L951" | variable == "p1" |
        variable == "padj_p1" ~ "Model 1",
      variable == "beta2" | variable == "U952"
      | variable == "L952" | variable == "p2" |
        variable == "padj_p2" ~ "Model 2",
      variable == "beta3" | variable == "U953"
      | variable == "L953" | variable == "p3" |
        variable == "padj_p3" ~ "Model 3")) %>% 
  # remove the numeric end of the colname (orig indicated model num)
  mutate(variable =  gsub("[1,2,3]", "", variable)) %>% 
  # widen data so cols for beta, L95, U95, p, padj
  pivot_wider(names_from = variable, values_from = value) %>% 
  mutate(Outcome = "Whole brain volume") %>%
  rename(coefficient = beta)

hv <- read.xlsx("results/metab/insight_metabolites.xlsx", sheet=2) %>%
  gather(variable, value, c("beta1", "beta2", "beta3", 
                            "L951", "L952", "L953",
                            "U951", "U952", "U953",
                            "p1", "p2", "p3",
                            "padj_p1", "padj_p2", "padj_p3")) %>%
  mutate(model = case_when(
    variable == "beta1" | variable ==  "U951"
    | variable == "L951" | variable == "p1" |
      variable == "padj_p1" ~ "Model 1",
    variable == "beta2" | variable == "U952"
    | variable == "L952" | variable == "p2" |
      variable == "padj_p2" ~ "Model 2",
    variable == "beta3" | variable == "U953"
    | variable == "L953" | variable == "p3" |
      variable == "padj_p3" ~ "Model 3")) %>% 
  mutate(variable =  gsub("[1,2,3]", "", variable)) %>% 
  pivot_wider(names_from = variable, values_from = value) %>% 
  mutate(Outcome = "Hippocampal volume") %>%
  rename(coefficient = beta)

am <- read.xlsx("results/metab/insight_metabolites.xlsx", sheet=3) %>%
  gather(variable, value, c("OR1", "OR2", "OR3", 
                            "L951", "L952", "L953",
                            "U951", "U952", "U953",
                            "p1", "p2", "p3",
                            "padj_p1", "padj_p2", "padj_p3")) %>%
  mutate(model = case_when(
    variable == "OR1" | variable ==  "U951"
    | variable == "L951" | variable == "p1" |
      variable == "padj_p1" ~ "Model 1",
    variable == "OR2" | variable == "U952"
    | variable == "L952" | variable == "p2" |
      variable == "padj_p2" ~ "Model 2",
    variable == "OR3" | variable == "U953"
    | variable == "L953" | variable == "p3" |
      variable == "padj_p3" ~ "Model 3")) %>% 
  mutate(variable =  gsub("[1,2,3]", "", variable)) %>% 
  pivot_wider(names_from = variable, values_from = value) %>% 
  mutate(Outcome = "Amyloid-beta status") %>%
  rename(coefficient = OR)

### Combine ----
fig3_hub <- rbind(bv, hv, am) %>%
  # filter for models 1 & 3
  filter(model == "Model 1" | model == "Model 3") %>%
  # indicate N and analysis type
  mutate(Analysis = "Imaging - Insight 46 (max N=437)", 
         # define effect direction & significance threshold
         Direction = factor(case_when(Outcome != "Amyloid-beta status" & coefficient <0  ~ "Negative",
                                      Outcome != "Amyloid-beta status" & coefficient >0  ~ "Positive",
                                      Outcome == "Amyloid-beta status" & coefficient <1 ~ "Positive", 
                                      Outcome == "Amyloid-beta status" & coefficient >1  ~ "Negative")),
         Significance = factor(case_when(
           padj_p <0.05 ~ "Adjusted threshold",
           padj_p >0.05 & p <0.05 ~ "Nominal",
           p >0.05  ~ "Non-sig"),
           levels = c("Adjusted threshold", "Nominal","Non-sig")),
         # text matrix with effect sizes
         textMatrix = ifelse(p<0.05, signif(coefficient,2), ""),
         # rename for plotting purposes
         model = case_when(model == "Model 1" ~ "Basic model",
                           model == "Model 3" ~ "Final model"),
         Outcome = factor(Outcome, levels = c("Whole brain volume", 
                                              "Hippocampal volume",
                                              "Amyloid-beta status"))) %>%
  # select cols of interest 
  select("BIOCHEMICAL", "SUPER.PATHWAY", "SUB.PATHWAY", 
         "module", "model", "p", "coefficient", "L95", "U95", 
         "Outcome", "Direction", "Significance", "Analysis", 
         "textMatrix")

### Format PGS results for plotting ----
fig3_pgs <- read.csv("results/pgs/PGS_metabs.csv") %>%
  select(BIOCHEMICAL,
    beta_APOE,beta_noAPOE,beta_APOE0.1,beta_noAPOE0.1,
    p_APOE,p_noAPOE,p_APOE0.1,p_noAPOE0.1,
    padj_APOE,padj_noAPOE,padj_APOE0.1,padj_noAPOE0.1) %>%
  # stack results
  gather(variable,value,
    c("beta_APOE","beta_noAPOE",
      "beta_APOE0.1","beta_noAPOE0.1",
      "p_APOE","p_noAPOE",
      "p_APOE0.1","p_noAPOE0.1",
      "padj_APOE","padj_noAPOE",
      "padj_APOE0.1","padj_noAPOE0.1")) %>% 
  # create col indicating model
  mutate(model = case_when(
    variable == "beta_APOE" | variable ==  "p_APOE" | 
      variable == "padj_APOE" | variable == "beta_APOE0.1" |
      variable ==  "p_APOE0.1" | variable == "padj_APOE0.1" ~ "APOE PRS",
    variable == "beta_noAPOE" | variable == "p_noAPOE" | 
      variable == "padj_noAPOE" | variable == "beta_noAPOE0.1" |
      variable == "p_noAPOE0.1" | variable == "padj_noAPOE0.1" ~ "No APOE PRS"),
    variable = gsub("_APOE", "", variable),
    variable = gsub("_noAPOE", "", variable)) %>%
  # widen data so cols for beta, L95, U95, p, padj
  pivot_wider(names_from = variable, values_from = value) 

# Select sig prs metabs
prs_mets <- fig3_pgs %>% 
  filter(p<0.05 | p0.1<0.05) %>% 
  pull(BIOCHEMICAL) %>% 
  as.character()

fig3_pgs %<>%
  # filter res to best threshold
  mutate(best_thresh = case_when(
    BIOCHEMICAL %in% prs_mets & p<p0.1 ~ "GW", p<0.05 |
      BIOCHEMICAL %in% prs_mets & p>p0.1 ~ "0.1"),
    beta_final = case_when(
      best_thresh == "0.1" & model == "APOE PRS"~ beta0.1,
      best_thresh == "0.1" & model == "No APOE PRS"~ beta0.1,
      best_thresh == "GW" & model == "APOE PRS"~ beta,
      best_thresh == "GW" & model == "No APOE PRS"~ beta,
      is.na(best_thresh) & model == "APOE PRS" ~ beta,
      is.na(best_thresh) & model == "No APOE PRS" ~ beta),
    p_final = case_when(
      best_thresh == "0.1" & model == "APOE PRS"~ p0.1,
      best_thresh == "0.1" & model == "No APOE PRS"~ p0.1,
      best_thresh == "GW" & model == "APOE PRS"~ p,
      best_thresh == "GW" & model == "No APOE PRS"~ p,
      is.na(best_thresh) & model == "APOE PRS" ~ p,
      is.na(best_thresh) & model == "No APOE PRS" ~ p),
    padj_final = case_when(
      best_thresh == "0.1" & model == "APOE PRS"~ padj0.1,
      best_thresh == "0.1" & model == "No APOE PRS"~ padj0.1,
      best_thresh == "GW" & model == "APOE PRS"~ padj,
      best_thresh == "GW" & model == "No APOE PRS"~ padj,
      is.na(best_thresh) & model == "APOE PRS" ~ padj,
      is.na(best_thresh) & model == "No APOE PRS" ~ padj)) %>%
  # get pathways
  left_join(pathway, "BIOCHEMICAL") %>%
  # get modules
  left_join(metab_MM, "BIOCHEMICAL") %>%
  # set significance threshold, effect direction, & coefficient text matrix
  mutate(Analysis=rep("NSHD (N=1638)"),
         Outcome = model,
         model = "AD PRS",
         L95 = "",
         U95 = "",
         Significance = case_when(
           padj_final<0.05 ~ "Adjusted threshold",
           padj_final >0.05 & p_final<0.05 ~ "Nominal",
           p_final>0.05 ~"Non-sig"),
         Direction = case_when(
           beta_final <0 ~ "Positive",
           beta_final>0 ~ "Negative"),
         textMatrix = ifelse(p_final<0.05, signif(beta_final,2), "")) %>%
  # select cols of interest
  select(BIOCHEMICAL, SUPER.PATHWAY, SUB.PATHWAY, 
         module, model, p=p_final, 
         coefficient=beta_final, L95, U95, 
         Outcome, Direction, Significance, 
         Analysis, textMatrix) 

## Combine data ----
fig3 <- rbind(fig3_hub %>% filter(BIOCHEMICAL %in% fig3_pgs$BIOCHEMICAL), # filter for 30 hubs followed up with pgs
              fig3_pgs) %>%
  # bold markdown for metabolites assoc in previous cognition analyses
  mutate(BIOCHEMICAL = case_when(
    BIOCHEMICAL %in% cog_hubs ~ paste0("**", BIOCHEMICAL, "**"), 
    TRUE ~ paste0(BIOCHEMICAL)),
    # set order for plotting
    module = factor(module, levels = c("blue", "brown", "yellow", 
                                       "cyan", "turquoise", "purple",
                                       "greenyellow", "magenta", "red", "tan")),
    Analysis = factor(Analysis, 
                      levels = c("Imaging - Insight 46 (max N=437)", "NSHD (N=1638)")),
    SUB.PATHWAY = case_when(
      BIOCHEMICAL == "**X - 21470**" ~ "Unknown",
      BIOCHEMICAL != "**X - 21470**" ~ SUB.PATHWAY),
    BIOCHEMICAL = case_when(BIOCHEMICAL == "**2,3-dihydroxy-5-methylthio-4-pentenoate (DMTPA)***" ~ "**DMTPA***",
                            BIOCHEMICAL != "**2,3-dihydroxy-5-methylthio-4-pentenoate (DMTPA)***" ~ BIOCHEMICAL))

# Order pathways by module
pathways <- fig3 %>% 
  arrange(module)  %>% 
  select(SUB.PATHWAY) %>% 
  unique() %>%
  pull(SUB.PATHWAY)

fig3 %<>% 
  mutate(SUB.PATHWAY = factor(SUB.PATHWAY, levels = pathways))

## Plot ----
fig3_plot = fig3 %>% 
  ggplot(aes(Outcome, BIOCHEMICAL, SUB.PATHWAY)) + 
  # heatmap tiles, fill = effect direction and transparency = p threshold
  geom_tile(aes(fill = Direction, alpha = Significance), color = NA) + 
  scale_fill_manual(values = c("red", "blue", "white"))+ 
  scale_alpha_manual(values = c(1, 0.5, 0)) + #setting transparency such that non-sig 0.05, nominal 0.5, bonf 1
  # no legends
  guides(fill = "none", alpha = "none") +
  scale_y_discrete() +
  # coefficients in tiles where p<0.05
  geom_text(aes(label = textMatrix), size=4.8, fontface="bold") + 
  # split by pathway, analysis and model model
  facet_nested(SUB.PATHWAY~ Analysis + model, 
               scales = "free", space = "free", switch = "x", nest_line=T) + 
  # remove axis titles
  ylab("") + xlab("") +
  # adjust text (render some as markdown) and remove background elements
  theme(panel.background = element_blank(), 
        panel.spacing.y = unit(0.04, 'lines'),
        strip.background.y  = element_rect(fill = NA),
        strip.text.y = element_text(angle=0, hjust = 0, margin = margin(l=18), size=15, face="bold"), 
        strip.text.x = element_text(face="bold", size=14),
        axis.text.x = element_text(angle=330, size=14.5, hjust=0, face="bold"),
        axis.text.y = element_markdown(size=14.5), 
        axis.title.x = element_text(margin = margin(t=15, r=0, b=0, l=0), face="bold"),
        axis.title.y = element_text(margin = margin(t=0, r=15, b=0, l=0), face="bold"))

# Add colour panel at RHS
fig3_plot <- as_grob(fig3_plot)
stripr <- which(grepl('strip-r', fig3_plot$layout$name))

# Indicate colours (number is no of pathways in each one)
fills <- c(rep("blue", 5), 
           rep("brown", 3), 
           rep("yellow", 1), 
           rep("turquoise", 2), 
           rep("magenta", 1), 
           rep("red", 1))

k <- 1
for (i in stripr) {
  lp <- linesGrob(x=unit(c(0,0),"npc"), y=unit(c(0,1),"npc"), 
                  # change facet border to just one side, very thick, and coloured by module
                  gp=gpar(col=paste(fills[k]), lwd=40)) 
  fig3_plot$grobs[[i]]$grobs[[1]]$children[[1]] <- lp 
  k <- k+1
}

# Save
ggsave(fig3_plot ,height=13.5, width=18, units="in", filename="figures/fig3plot.pdf", dpi=500)
