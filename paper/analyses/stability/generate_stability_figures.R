library(readr)
library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)

## Load in Test Case Data

test_cases_df <- read_excel("../../../tests/test_cases.xlsx", sheet = "generated_seqs") %>% 
  select(antigen_id, model, seq_id)

## Load in FoldX Results and Join

foldx_df <- read_csv("foldx_results.csv") %>% 
  rename(
    total_energy = `Total Energy`,
    solvation_polar = `Solvation Polar`,
    solvation_hydrophobic = `Solvation Hydrophobic`
    ) %>% 
  select(seq_id, total_energy, solvation_polar,solvation_hydrophobic) %>%
  inner_join(test_cases_df, by = "seq_id")

## Total Energy
te_mu <- plyr::ddply(foldx_df, c("antigen_id", "model"), summarise, grp.mean=mean(total_energy))

te_plot <- ggplot(foldx_df, aes(x=total_energy, fill=model)) +
  geom_density(alpha=0.4) +
  geom_vline(data=te_mu, aes(xintercept=grp.mean, color=model),
             linetype="dashed") +
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) +
  facet_wrap(
    ~antigen_id,
    ncol = 7,
    scales = "free_x"
    ) +
  xlab("Total Energy (kcal/mol)") +
  theme_linedraw() +
  theme(
    legend.position="bottom",
    axis.text.x = element_text(angle = 90, hjust = 1)
  )
  
## Polar Solvation
ps_mu <- plyr::ddply(foldx_df, c("antigen_id", "model"), summarise, grp.mean=mean(solvation_polar))

ps_plot <- ggplot(foldx_df, aes(x=solvation_polar, fill=model)) +
  geom_density(alpha=0.4) +
  geom_vline(data=ps_mu, aes(xintercept=grp.mean, color=model),
             linetype="dashed") +
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) +
  facet_wrap(
    ~antigen_id,
    ncol = 7,
    scales = "free_x"
  ) +
  # xlim(-700, -500) + 
  xlab("Polar Solvation (kcal/mol)") +
  theme_linedraw() +
  theme(
    legend.position="bottom",
    axis.text.x = element_text(angle = 90, hjust = 1)
  )



## Hydrophobic Solvation
hs_mu <- plyr::ddply(foldx_df, c("antigen_id", "model"), summarise, grp.mean=mean(solvation_hydrophobic))

hs_plot <- ggplot(foldx_df, aes(x=solvation_hydrophobic, fill=model)) +
  geom_density(alpha=0.4) +
  geom_vline(data=hs_mu, aes(xintercept=grp.mean, color=model),
             linetype="dashed") +
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) +
  facet_wrap(
    ~antigen_id,
    ncol = 7,
    scales = "free_x"
  ) +
  xlab("Hydrophobic Solvation (kcal/mol)") +
  theme_linedraw() +
  theme(
    legend.position="bottom",
    axis.text.x = element_text(angle = 90, hjust = 1)
    )



stability_figure <- ggarrange(te_plot,
                              ps_plot,
                              hs_plot,
                              labels = c("A", "B", "C"),
                              ncol = 1, nrow = 3,
                              common.legend = TRUE,
                              legend = "bottom")

stability_figure
