library(readr)
library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)

## Remove outliers function
## Adapted from: https://stackoverflow.com/questions/49982794/remove-outliers-by-group-in-r
remove_outliers_iqr <- function(x) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = TRUE)
  H <- 1.5 * IQR(x, na.rm = TRUE)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  return(y)
}


## Load Data
test_cases_df <- read_excel("../tests/test_cases.xlsx", sheet = "generated_seqs")
affinity_df <- read_csv("../tests/binding_affinity_results.csv") %>% 
  mutate(
    seq_id = str_replace_all(seq_id, c("predicted__" = "", "__esm3" = "")),
    score = remove_outliers_iqr(score)
  ) %>% 
  filter(!seq_id %in% c("peleke-phi4_PD-1_0a", "peleke-phi4_BHRF1_0a", "peleke-phi4_PD-1_0a", "peleke-phi4_MBP_0a")) %>%
  filter(!is.na(score))


## Merge Data
affinity_merged_df <- test_cases_df %>%
  inner_join(affinity_df, by = "seq_id")


## Hits by antigen/model
threshold = -25
hits_by_antigen <- affinity_merged_df %>%
  group_by(antigen_id, model) %>%
  summarise(
    percentage_hits = mean(vdw < threshold) * 100,
    median_score = median(vdw)
    ) %>% 
  tidyr::pivot_wider(names_from = model, values_from = c(percentage_hits, median_score))


## Plotting

ggplot(data=affinity_merged_df, aes(x=score, group=model, fill=model)) +
  geom_density(adjust=1.5, alpha=0.6) +
  facet_wrap(~antigen_id) +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    axis.ticks.x=element_blank()
  ) +
  xlim(-250, 250)





ggdensity(affinity_merged_df, x = "score",
          add = "median", rug = TRUE,
          color = "model", fill = "model", facet.by = "antigen_id",
          palette = c("#00AFBB", "#E7B800", "#FC4E07")) + xlim(-250, 0)


## ggpubr boxplot

ggboxplot(affinity_merged_df,
          x = "model", y = "score",
          color = "model",
          palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          add = "jitter",
          facet.by = "antigen_id"
          ) + ylim(-250, 250)
