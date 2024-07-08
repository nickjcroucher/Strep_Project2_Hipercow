
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

# PLOT Stats_joined N50, colours based on delete_priority
col_map <- c("0" = "lightgreen", "1" = "red", "2" = "orange", "3" = "gold2")
stats_joined$col <- col_map[as.character(stats_joined$delete_priority)]

png("pictures/tree_pruning.png", width = 12, height = 18, unit = "cm", res = 1200)
par(mfrow = c(3,2), mar = c(3, 3, 2, 2), mgp = c(1.7, 0.7, 0), bty = "n")
plot(stats_joined$N50, stats_joined$L50,
     col = stats_joined$col,
     pch = 1, xlab = "", ylab = "L50")
plot(stats_joined$N50, stats_joined$mean,
     col = stats_joined$col,
     pch = 1, xlab = "", ylab = "mean")
plot(stats_joined$N50, stats_joined$median,
     col = stats_joined$col,
     pch = 1, xlab = "", ylab = "median")

plot(stats_joined$N50, stats_joined$sequence_count,
     col = stats_joined$col,
     pch = 1, xlab = "", ylab = "Sequence Count")
plot(stats_joined$N50, stats_joined$shortest,
     col = stats_joined$col,
     pch = 1, xlab = "", ylab = "Shortest")
plot(stats_joined$N50, stats_joined$total_bps,
     col = stats_joined$col,
     pch = 1, xlab = "", ylab = "Total BPS")

dev.off()

par(mfrow = c(1,1))
