
library(tidyverse)

stats <- read.csv("raw_data/gubbins/stats_compiled.csv")
stats$contigs <- substr(stats$file, 7, 51)

# Pruning based on tree output from microreact (https://microreact.org/upload)
suppressWarnings({
pr_1 <- read.table("raw_data/gubbins/leaf_labels_priority_1.txt", header = F)
pr_2 <- read.table("raw_data/gubbins/leaf_labels_priority_2.txt", header = F)
pr_3 <- read.table("raw_data/gubbins/leaf_labels_priority_3.txt", header = F)
})

pr_1$delete_priority <- 1
pr_2$delete_priority <- 2
pr_3$delete_priority <- 3

pr_all <- rbind(pr_1, pr_2, pr_3)

stats_joined <- dplyr::left_join(stats, pr_all, by = c("contigs" = "V1"))
stats_joined$delete_priority[is.na(stats_joined$delete_priority)] <- 0

stats_G <- stats_joined %>% 
  dplyr::group_by(delete_priority) %>% 
  dplyr::summarise(mean_L50 = mean(L50),
                   mean_N50 = mean(N50),
                   min_N50 = min(N50),
                   max_N50 = max(N50),
                   mean_gc_content = mean(gc_content),
                   mean_mean = mean(mean),
                   mean_seq_count = mean(sequence_count),
                   mean_total_bps = mean(total_bps)) %>% 
  dplyr::ungroup()

