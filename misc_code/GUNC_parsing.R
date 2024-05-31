################################################################################
### Code for parsing GUNC output
### Raphael Eisenhofer 2024
################################################################################

library(tidyverse)

df <- read_delim("data/GUNC.progenomes_2.1.maxCSS_level.tsv") %>%
  filter(!str_detect(genome, "ALL_dreped")) %>%
  mutate(overall_strategy = case_when(
    str_detect(genome, "individual") ~ "single_coverage",
    str_detect(genome, "cobinning_long") ~ "multicoverage_animal",
    str_detect(genome, "cobinning_treat") ~ "multicoverage_timepoint_all",
    str_detect(genome, "cobinning_cage_treat") ~ "multicoverage_timepoint_cage",
    str_detect(genome, "coassembly_long") ~ "coassembly_animal",
    str_detect(genome, "coassembly_treat") ~ "coassembly_timepoint_all",
    str_detect(genome, "coassembly_cage_treat") ~ "coassembly_timepoint_cage",
    str_detect(genome, "vamb_long") ~ "multisplit_animal",
    str_detect(genome, "vamb_treat") ~ "multisplit_timepoint_all",
    str_detect(genome, "vamb_cage_treat") ~ "multisplit_timepoint_cage",
    TRUE ~ genome))


#Plotting
#How many bins pass GUNC threshold per strategy?
df %>% 
  summarise(n = n(), .by = c(overall_strategy, pass.GUNC)) %>%
  ggplot(aes(x = overall_strategy, y = n, fill = pass.GUNC)) +
  geom_bar(stat = "identity", position = position_stack(reverse = FALSE)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, size = 14, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold")) +
  labs(y = "Number of bins", x = "Overall strategy")
  

#N contigs per strategy?
df %>% 
  filter(strategy != "cobinning_treat") %>%
  ggplot(aes(x = strategy, y = n_contigs)) +
  geom_boxplot() + 
  geom_jitter(width = 0.2, height = 0, alpha = 0.5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 14, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold")) +
  labs(y = "n_contigs per bin")



#Save summary for another figure

write_delim(df %>% select(genome, n_contigs, pass.GUNC, overall_strategy), "data/gunc_summary.tsv")
