################################################################################
### Code for parsing GUNC output
### Raphael Eisenhofer 2024
################################################################################

library(tidyverse)

df <- read_delim("data/GUNC.progenomes_2.1.maxCSS_level.tsv") %>%
  mutate(strategy = case_when(str_detect(genome, "ALL_dreped") ~ "all_dereplicated",
                              str_detect(genome, "coassembly_cage_treat") ~ "coassembly_cage_treat",
                              str_detect(genome, "coassembly_treat") ~ "coassembly_treat",
                              str_detect(genome, "coassembly_long") ~ "coassembly_long",
                              str_detect(genome, "cobinning_cage_treat") ~ "cobinning_cage_treat",
                              str_detect(genome, "cobinning_treat") ~ "cobinning_treat",
                              str_detect(genome, "cobinning_long") ~ "cobinning_long",
                              str_detect(genome, "vamb_cage_treat") ~ "vamb_cage_treat",
                              str_detect(genome, "vamb_treat") ~ "vamb_treat",
                              str_detect(genome, "vamb_long") ~ "vamb_long",
                              str_detect(genome, "individual") ~ "individual",
                              ),
         )


#Plotting
#How many bins pass GUNC threshold per strategy?
df %>% 
  filter(strategy != "cobinning_treat") %>%
  summarise(n = n(), .by = c(strategy, pass.GUNC)) %>%
  ggplot(aes(x = strategy, y = n, fill = pass.GUNC)) +
  geom_bar(stat = "identity", position = position_stack(reverse = FALSE)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 14, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold")) +
  labs(y = "n bins")
  

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

write_delim(df %>% select(genome, n_contigs, pass.GUNC), "data/gunc_summary.tsv")
