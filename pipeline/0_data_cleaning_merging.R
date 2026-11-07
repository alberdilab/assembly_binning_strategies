### Cleaning/pivotting data
### Raphael Eisenhofer 2023

library(tidyverse)
library(janitor)


# Consolidate assembly statistics

coassembly_cage_treat_files <- list.files("data/tables/", 
                                          pattern = "coassembly_cage_treat",
                                          full.names = T)
coassembly_long_files <- list.files("data/tables/", 
                                          pattern = "coassembly_long",
                                          full.names = T)
coassembly_treat_files <- list.files("data/tables/", 
                                          pattern = "coassembly_treat",
                                          full.names = T)

coassembly_cage_treat <- map(coassembly_cage_treat_files, ~ read_delim(file = .x)) %>%
  list_rbind() %>%
  mutate(assembly = "coassembly",
         strategy = "cage_treatment",
         sample = str_replace_all(sample, "_subsam", ""))

coassembly_long <- map(coassembly_long_files, ~ read_delim(file = .x)) %>%
  list_rbind() %>%
  mutate(assembly = "coassembly",
         strategy = "longitudinal",
         sample = str_replace_all(sample, "_subsam", ""))

coassembly_treat <- map(coassembly_treat_files, ~ read_delim(file = .x)) %>%
  list_rbind() %>%
  mutate(assembly = "coassembly",
         strategy = "treatment",
         sample = str_replace_all(sample, "_subsam", ""))

individual_assembly <- read_delim("data/tables/assembly_summary.tsv.gz") %>%
  mutate(assembly = "individual",
         strategy = "individual",
         sample = str_replace_all(sample, "_subsam", "")) %>%
  rename(assembly_mapping_percent = aseembly_mapping_percent)

df_list <- list(coassembly_cage_treat, coassembly_long, coassembly_treat, individual_assembly)

assembly_comparison <- df_list %>% reduce(full_join)

assembly_comparison %>%
  ggplot(aes(y = assembly_mapping_percent, x = sample)) +
  geom_point() +
  facet_grid(~strategy)

assembly_comparison %>%
  ggplot(aes(y = N50, x = sample)) +
  geom_point() +
  facet_grid(~strategy)

assembly_comparison %>%
  ggplot(aes(x = .panel_x, y = .panel_y)) +
  geom_point(aes(color = assembly), alpha = .2) +
  geom_smooth(aes(color = assembly), method = "lm") +
  ggforce::geom_autodensity(aes(color = assembly, fill = after_scale(color)), alpha = .7) +
  scale_color_brewer(palette = "Set2", name = NULL) +
  ggforce::facet_matrix(vars(N50, L50, num_contigs, largest_contig, total_length, assembly_mapping_percent), layer.lower = 2, layer.diag = 3)


assembly_comparison %>%
  filter(assembly == "coassembly") %>%
  ggplot(aes(x = .panel_x, y = .panel_y)) +
  geom_point(aes(color = strategy), alpha = .2) +
  geom_smooth(aes(color = strategy), method = "lm") +
  ggforce::geom_autodensity(aes(color = strategy, fill = after_scale(color)), alpha = .7) +
  scale_color_brewer(palette = "Set2", name = NULL) +
  ggforce::facet_matrix(vars(N50, L50, num_contigs, largest_contig, total_length, assembly_mapping_percent), layer.lower = 2, layer.diag = 3)


######################################
# Consolidate binning results
######################################

## declare strategies
treatments <- read_delim("data/treatments.tsv", delim = '\t', col_names = "treatment")
cage_treatments <- read_delim("data/cage_treat_folders.tsv", delim = '\t', col_names = "cage_treatment")
mice <- read_delim("data/mice.tsv", delim = '\t', col_names = "mouse")

## import n_contigs
ncontigs <- read_delim("data/tables/ncontigs.tsv", delim = '\t', col_names = c("bin_id", "n_contigs")) %>%
  mutate(bin_id = str_replace_all(bin_id, ".fa.gz", ""))

######################################
## Individual assembly/binning MAGs
######################################
df_individual <- read_delim("data/tables/individual_metawrap_70_10_bins.stats.gz") %>%
  mutate(assembly = "individual", strategy = "individual_binning",
         genome = str_replace_all(genome, "^", "individual_")) %>%
  rename(bin = genome)


######################################
## coassembly treatment MAGs
######################################

#create function for reading input files
read_data <- function(x) {
  read_delim(file = paste0('data/tables/', x, '_metawrap_70_10_bins.stats.gz'), delim = '\t') %>%
    mutate(assembly = "coassembly", strategy = "treatment",
           bin = str_replace_all(bin, "^", "coassembly_treat_"))
}

#iterate using purrr's map(), then merge into a single dataframe
df_coassembly_treatment <- map(treatments$treatment, read_data) %>%
  list_rbind()


######################################
## coassembly cage_treatment MAGs
######################################

#create function for reading input files
read_data1 <- function(x) {
  read_delim(file = paste0('data/tables/', x, '_metawrap_70_10_bins.stats.gz'), delim = '\t') %>%
    mutate(assembly = "coassembly", strategy = "cage_treatment",
           bin = str_replace_all(bin, "^", "coassembly_cage_treat_"))
}

#iterate using purrr's map(), then merge into a single dataframe
df_coassembly_cage_treatment <- map(cage_treatments$cage_treatment, read_data1) %>%
  list_rbind()


######################################
## coassembly longitudinal MAGs
######################################

#create function for reading input files
read_data2 <- function(x) {
  read_delim(file = paste0('data/tables/', x, '_metawrap_70_10_bins.stats.gz'), delim = '\t') %>%
    mutate(assembly = "coassembly", strategy = "longitudinal",
           bin = str_replace_all(bin, "^", "coassembly_long_"))
}

#iterate using purrr's map(), then merge into a single dataframe
df_coassembly_longitudinal <- map(mice$mouse, read_data2) %>%
  list_rbind()



######################################
## cobinning treatment MAGs
######################################
cobinning_treat_files <- list.files("data/tables/", 
                                     pattern = "cobinning_treat",
                                     full.names = T)

df_cobinning_treatment <- map(cobinning_treat_files, ~ read_delim(file = .x)) %>%
  list_rbind() %>%
  mutate(assembly = "cobinning",
         strategy = "treatment",
         bin = str_replace_all(bin, "^", "cobinning_treat_"))


######################################
## cobinning cage_treatment MAGs
######################################
cobinning_cage_treatment_files <- list.files("data/tables/", 
                                    pattern = "cobinning_cage_treat",
                                    full.names = T)

df_cobinning_cage_treatment <- map(cobinning_cage_treatment_files, ~ read_delim(file = .x)) %>%
  list_rbind() %>%
  mutate(assembly = "cobinning",
         strategy = "cage_treatment",
         bin = str_replace_all(bin, "^", "cobinning_cage_treat_"))


######################################
## cobinning longitudinal MAGs
######################################
cobinning_long_files <- list.files("data/tables/", 
                                    pattern = "cobinning_long",
                                    full.names = T)

df_cobinning_long <- map(cobinning_long_files, ~ read_delim(file = .x)) %>%
  list_rbind() %>%
  mutate(assembly = "cobinning",
         strategy = "longitudinal",
         bin = str_replace_all(bin, "^", "cobinning_long_"))


######################################
## vamb treatment MAGs
######################################

# Load in vamb mag stats
vamb_stats <- read_delim("data/tables/vamb_mag_stats.tsv.gz") %>%
  mutate(filename = str_replace_all(filename, "/maps/.*vamb_long.*vamb_bins/", "vamb_long_"),
         filename = str_replace_all(filename, "/maps/.*vamb_cage_treatment.*vamb_bins/", "vamb_cage_treat_"),
         filename = str_replace_all(filename, "/maps/.*vamb_treat.*vamb_bins/", "vamb_treat_"),
         bin_id = str_replace_all(filename, ".fa", "")) %>%
  select(bin_id, n_contigs, scaf_L50, scaf_bp, gc_avg) %>%
  rename(size = scaf_bp, GC = gc_avg, N50 = scaf_L50)


vamb_treatment_files <- list.files("data/tables/", 
                                   pattern = "vamb_treatment",
                                   full.names = T)

df_vamb_treatment <- map(vamb_treatment_files, ~ read_delim(file = .x)) %>%
  list_rbind() %>%
  clean_names() %>%
  mutate(assembly = "multi_split",
         strategy = "treatment",
         bin_id = str_replace_all(bin_id, "^", "vamb_treat_"),
         bin_id = str_replace_all(bin_id, "_OP_OP", "_OP")) %>%
  select(bin_id, completeness, contamination, assembly, strategy) %>%
  left_join(., vamb_stats, by = join_by(bin_id)) %>%
  filter(completeness >= 70 & contamination <= 10)


######################################
## vamb cage_treatment MAGs
######################################
vamb_cage_treatment_files <- list.files("data/tables/", 
                                   pattern = "vamb_cage_treatment",
                                   full.names = T)

df_vamb_cage_treatment <- map(vamb_cage_treatment_files, ~ read_delim(file = .x)) %>%
  list_rbind() %>%
  clean_names() %>%
  mutate(assembly = "multi_split",
         strategy = "cage_treatment",
         bin_id = str_replace_all(bin_id, "^", "vamb_cage_treat_")) %>%
#         bin_id = str_replace_all(bin_id, "_C.*_C", "_C")) %>%
  select(bin_id, completeness, contamination, assembly, strategy) %>%
  left_join(., vamb_stats, by = join_by(bin_id)) %>%
  filter(completeness >= 70 & contamination <= 10)


######################################
## vamb longitudinal MAGs
######################################
vamb_longitudinal_files <- list.files("data/tables/", 
                                   pattern = "vamb_long",
                                   full.names = T)

df_vamb_longitudinal <- map(vamb_longitudinal_files, ~ read_delim(file = .x)) %>%
  list_rbind() %>%
  clean_names() %>%
  mutate(assembly = "multi_split",
         strategy = "longitudinal",
         bin_id = str_replace_all(bin_id, "^", "vamb_long_"),
         bin_id = str_replace_all(bin_id, "_C.*_C", "_C")) %>%
  select(bin_id, completeness, contamination, assembly, strategy) %>%
  left_join(., vamb_stats, by = join_by(bin_id)) %>%
  filter(completeness >= 70 & contamination <= 10)


######################################
## Combine into one dataframe
######################################

#load in derep clustering info
drep_clustering <- read_delim("data/drep_all/Cdb.csv.gz") %>%
  mutate(genome = str_replace_all(genome, ".fa.gz$", ""),
         genome = str_replace_all(genome, ".fa$", ""),
         genome = str_replace_all(genome, "C11F2_C11F2", "C11F2"),
         genome = str_replace_all(genome, "C06M3_C06M3", "C06M3"),
         genome = str_replace_all(genome, "C06M5_C06M5", "C06M5"),
         genome = str_replace_all(genome, "OP_OP", "OP")) %>%
  rename(bin_id = genome) %>%
  select(bin_id, secondary_cluster, primary_cluster)


all_df <- rbind(df_coassembly_cage_treatment, df_coassembly_longitudinal, df_coassembly_treatment,
      df_cobinning_cage_treatment, df_cobinning_long, df_cobinning_treatment,
      df_individual) %>%
  rename(bin_id = bin) %>%
  select(!binner) %>%
  select(!lineage) %>%
  inner_join(., ncontigs, by = join_by(bin_id)) %>%
  rbind(., df_vamb_cage_treatment, df_vamb_longitudinal, df_vamb_treatment) %>%
  inner_join(., drep_clustering, by = join_by(bin_id)) %>%
  unite(col = "overall_strategy", assembly, strategy, sep = "_", remove = F)


write_delim(all_df, "data/binning_stats.tsv")




