#====================#
#      Figure 2      #
#====================#

## Set up ----
filepath <- ""
setwd(filepath)

### Load libraries ----
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(gtable)
library(ggtext)
library(magrittr)
library(openxlsx)
library(tidyr)

## Figure 2a (module heatmap) ----
### Prep data ----
#### Brain volume ----
bv <- read.xlsx("results/module/insight_modules.xlsx", sheet = 1) %>%
  mutate(outcome = "Whole brain volume") %>%
  # stack results
  gather(variable, value,
    c("beta1","beta2","beta3",
      "L951","L952","L953",
      "U951","U952","U953",
      "p1","p2","p3",
      "padj_p1","padj_p2","padj_p3")) %>%
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
        variable == "padj_p3" ~ "Model 3"),
    # remove the numeric end of the colname (orig indicated model num)
    variable =  gsub("[1,2,3]", "", variable)) %>%
  # widen data so cols for beta, L95, U95, p, padj
  pivot_wider(names_from = variable, values_from = value)

#### Hippocampal volume ----
hv <- read.xlsx("results/module/insight_modules.xlsx", sheet = 2) %>%
  mutate(outcome = "Hippocampal volume") %>%
  gather(variable, value,
         c("beta1","beta2","beta3",
           "L951","L952","L953",
           "U951","U952","U953",
           "p1","p2","p3",
           "padj_p1","padj_p2","padj_p3")) %>%
  mutate(model = case_when(
    variable == "beta1" | variable ==  "U951"
    | variable == "L951" | variable == "p1" |
      variable == "padj_p1" ~ "Model 1",
    variable == "beta2" | variable == "U952"
    | variable == "L952" | variable == "p2" |
      variable == "padj_p2" ~ "Model 2",
    variable == "beta3" | variable == "U953"
    | variable == "L953" | variable == "p3" |
      variable == "padj_p3" ~ "Model 3"),
    variable =  gsub("[1,2,3]", "", variable)) %>%
  pivot_wider(names_from = variable, values_from = value)

#### Amyloid status ----
am <- read.xlsx("results/module/insight_modules.xlsx", sheet = 3) %>%
  mutate(outcome = paste(expression("Amyloid-\u03B2 status"))) %>%
  gather(variable, value,
         c("OR1","OR2","OR3",
           "L951","L952","L953",
           "U951","U952","U953",
           "p1","p2","p3",
           "padj_p1","padj_p2","padj_p3")) %>%
  # create col indicating model
  mutate(model = case_when(
    variable == "OR1" | variable ==  "U951"
    | variable == "L951" | variable == "p1" |
      variable == "padj_p1" ~ "Model 1",
    variable == "OR2" | variable == "U952"
    | variable == "L952" | variable == "p2" |
      variable == "padj_p2" ~ "Model 2",
    variable == "OR3" | variable == "U953"
    | variable == "L953" | variable == "p3" |
      variable == "padj_p3" ~ "Model 3"),
    # remove the numeric end of the colname (orig indicated model num)
    variable =  gsub("[1,2,3]", "", variable)) %>%
  # widen data so cols for beta, L95, U95, p, padj
  pivot_wider(names_from = variable, values_from = value) %>%
  # Convert OR to effect direction
  mutate(Direction = factor(case_when(
    OR > 1 & p < 0.05 ~ "Negative",
    OR < 1 & p < 0.05 ~ "Positive",
    p > 0.05 ~ "NA"),
    levels = c("Negative", "Positive", "NA"))) %>%
  # rename to coefficient so can combine with other res (with betas)
  rename("coefficient" = "OR")

#### Combine data ----
# Order of modules for plotting 
# (ordered based on split in dendrogram)
module <- c("cyan","green","black",
            "brown","greenyellow","yellow",
            "pink","red","turquoise",
            "tan","salmon","blue",
            "magenta","purple")

fig2a <- rbind(bv, hv) %>%
  # convert betas to effect direction
  mutate(Direction = factor(case_when(
      beta < 0 & p < 0.05 ~ "Negative",
      beta > 0 & p < 0.05 ~ "Positive",
      p > 0.05 ~ "NA"), 
    levels = c("Negative", "Positive", "NA"))) %>%
  # rename to coefficient so can combine with amyloid res
  rename("coefficient" = "beta")  %>%
  # combine with amyloid
  rbind(am) %>%
  # define significance threshold & set module order for plotting
  mutate(threshold = factor(case_when(
        padj_p < 0.05 ~ "Adjusted threshold",
        padj_p > 0.05 & p < 0.05 ~ "Nominal",
        p > 0.05 ~ "Non-threshold"),
      levels = c("Adjusted threshold", "Nominal", "Non-threshold")),
    Module = factor(Module, levels = module),
    model = as.character(model)) %>%
  # filter for models 1 & 3
  filter(model == "Model 1" | model == "Model 3") %>%
  # rename for plotting purposes
  mutate(model = case_when(
    model == "Model 1" ~ "Basic model",
    model == "Model 3" ~ "Final model"),
    # text matrix with effect sizes
    textMatrix = ifelse(p < 0.05, paste(signif(coefficient, 2)), ""))

# Set coordinates for module coloured rectangles at top of plot
fig2a %<>% mutate(
    xmin = case_when(
      Module == "cyan" ~ 0.5,
      Module == "green" ~ 1.5,
      Module == "black" ~ 2.5,
      Module == "brown" ~ 3.5,
      Module == "greenyellow" ~ 4.5,
      Module == "yellow" ~ 5.5,
      Module == "pink" ~ 6.5,
      Module == "red" ~ 7.5,
      Module == "turquoise" ~ 8.5,
      Module == "tan" ~ 9.5,
      Module == "salmon" ~ 10.5,
      Module == "blue" ~ 11.5,
      Module == "magenta" ~ 12.5,
      Module == "purple" ~ 13.5),
    xmax = case_when(
      Module == "cyan" ~ 1.5,
      Module == "green" ~ 2.5,
      Module == "black" ~ 3.5,
      Module == "brown" ~ 4.5,
      Module == "greenyellow" ~ 5.5,
      Module == "yellow" ~ 6.5,
      Module == "pink" ~ 7.5,
      Module == "red" ~ 8.5,
      Module == "turquoise" ~ 9.5,
      Module == "tan" ~ 10.5,
      Module == "salmon" ~ 11.5,
      Module == "blue" ~ 12.5,
      Module == "magenta" ~ 13.5,
      Module == "purple" ~ 14.5
    ),
    ymin = 3.6, ymax = 4.1)

### Plot ----
fig2aplot = fig2a %>%
  ggplot(aes(Module, outcome)) +
  # heatmap tiles, fill = effect direction and transparency = p threshold
  geom_tile(aes(fill = Direction, alpha = threshold), color = NA) +
  scale_fill_manual(values = c("red", "blue", "white")) +
  scale_alpha_manual(values = c(1, 0.5, 0)) +
  # no legends
  guides(fill = "none", alpha = "none") +
  # coefficients in tiles where p<0.05
  geom_text(aes(label = textMatrix), size = 5.5) +
  # module names with N metabs in module 
  # + bold markdown for those assoc in previous cognition analyses
  scale_x_discrete(
    position = "top",
    limits = module,
    breaks = module,
    labels = c(
      "**Cyan<br>(n=22)**",
      "Green<br>(n=59)",
      "Black<br>(n=45)",
      "**Brown<br>(n=130)**",
      "Green-<br>yellow<br>(n=30)",
      "**Yellow<br>(n=70)**",
      "Pink<br>(n=35)",
      "Red<br>(n=47)",
      "**Turquoise<br>(n=192)**",
      "Tan<br>(n=26)",
      "Salmon<br>(n=25)",
      "Blue<br>(n=132)",
      "Magenta<br>(n=34)",
      "**Purple<br>(n=33)**")) +
  # module coloured rectangles at top of plot
  geom_rect(aes(xmin = xmin, xmax = xmax,
    ymin = ymin, ymax = ymax), fill = fig2a$Module) +
    # remove axis titles
  ylab("") + xlab("") +
  # split by model
  facet_grid(model ~ .,
             scales = "free_y",
             space = "free_y",
             switch = "y") +
  # adjust text (render some as markdown) and remove background elements
    theme(strip.text.x = element_blank(),
          strip.text.y = element_text(size = 16, face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(face = "bold"),
          axis.text.x.top = element_markdown(size = 13.5),
          axis.text.y = element_text(size = 15),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_rect(fill = "white"))

# Save
ggsave(plot = fig2aplot,
      filename = "figures/fig2aplot.png", 
      height = 5,
      width = 13,
      dpi = 500)

## Figure 2b (module pathways) ----
### Read in ORA pathway data ----
pathway <- readRDS("dat/other/ORA_pathway.rds") %>%
  # filter for significant pathways
  filter(padj < 0.05)

### Create dataframe for plotting ----
# No pathways in tan module, so restrict to non tan modules
ORA_mod <- module[module != "tan"]

# Create df for the module pathways
fig2b <- data.frame(row.names = ORA_mod, 
                    Pathway = rep(NA, length(ORA_mod)))

# Concat pathways for each module
for (i in ORA_mod){
  fig2b[i, 1] <- paste(pathway %>%  
                         filter(module == i) %>% 
                         arrange(padj) %>% 
                         pull(pathway) %>% 
                         paste(collapse='; '))
}

fig2b %<>%
  # add module column
  mutate(Module = row.names(fig2b)) %>%
  select(Module, Pathway) %>%
  # add in line breaks for long lines
  mutate(Pathway = case_when(Module == "turquoise" ~ "Histidine Metabolism; Urea cycle; Arginine and Proline Metabolism;\nMethionine, Cysteine, SAM and Taurine Metabolism; Tyrosine Metabolism; Polyamine Metabolism;\n Gamma-glutamyl Amino Acid; Leucine, Isoleucine and Valine Metabolism; Lysine Metabolism",
                             Module == "blue" ~ "Long Chain Polyunsaturated Fatty Acid (n3 and n6); Long Chain Fatty Acids;\nMedium Chain Fatty Acid; Fatty Acid, Monohydroxy; Fatty Acid, Dicarboxylate;\nEndocannabinoid; Fatty Acid, Branched",
                             TRUE ~ Pathway))

### Plot ----
fig2bplot <- tableGrob(fig2b, rows = NULL) %>%
  gtable_add_grob(grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                     t = 2, b = nrow(.), l = 1, r = ncol(.)) %>%
  gtable_add_grob(grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                     t = 1, l = 1, r = ncol(.))

grid.draw(fig2bplot)

# Save
ggsave(plot=fig2bplot, filename="figures/fig2bplot.pdf", h=10, w=10)

# End ----

