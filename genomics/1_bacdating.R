# devtools::install_github("xavierdidelot/BactDating")

library(tidyverse)
library(BactDating)
library(ape)
library(readxl)
library(coda)
library(socialmixr)

source("global/all_function.R") # Collected functions stored here!

# see help(package='BactDating') for more info
# Time-scaled tree with BactDating
# https://xavierdidelot.github.io/BactDating/articles/Staph.html

run_bacdating <- function(nbIts){
  # 1. Data wrangling ##########################################################
  # library(data.table)
  # data <- fread("raw_data/gubbins/spn_uk_dude.csv")
  data <- readxl::read_excel("raw_data/gubbins/ukhsa_assemblies_02_07_24.xlsx")
  dat <- readxl::read_excel("raw_data/serotype1_UKHSA_imperial_date_age_region_MOLIS_withdeath_meningitis_clean.xlsx") #%>% 
  
  tre <- BactDating::loadGubbins("raw_data/gubbins/n703/n703_")
  
  
  ##############################################################################
  tre_names <- as.data.frame(tre$tip.label) #%>% 
  # rename(ID = 'tre$tip.label')
  tre_names$ID <- substr(tre$tip.label, 1, 8)
  tre_names <- dplyr::left_join(tre_names, data, by = c("ID" = "ngsid"))
  tre_names <- dplyr::left_join(tre_names, dat, by = c("ID.y" = "ID"))
  tre_names <- tre_names %>% 
    dplyr::mutate(Earliestspecimendate = as.Date(Earliestspecimendate),
                  year = year(Earliestspecimendate))
  
  # 2. BacDating ###############################################################
  # https://xavierdidelot.github.io/BactDating/articles/yourData.html
  d <- cbind(tre_names$year, tre_names$year+1) # 2-d matrix according to the articles above
  
  set.seed(0)
  
  dir.create("outputs/genomics/choosen_n703", FALSE, TRUE)
  dir.create("pictures/genomics/choosen_n703", FALSE, TRUE)
  
  res_pr <- BactDating::bactdate(tre,d,nbIts=nbIts, # Put 1e6 or 1e10 on hipercow
                                 model = "strictgamma", # "arc",
                                 showProgress = T)
  
  saveRDS(res_pr, "outputs/genomics/choosen_n703/mcmc_bacdating.rds")
  
  # Figures!
  png("pictures/genomics/choosen_n703/tree_treeCI.png", width = 24, height = 12, unit = "cm", res = 1200)
  par(mfrow = c(1,1), mar = c(3, 3, 2, 2), mgp = c(1.7, 0.7, 0), bty = "n")
  plot(res_pr,'treeCI',show.tip.label = F)
  dev.off()
  
  png("pictures/genomics/choosen_n703/tree_trace1.png", width = 24, height = 12, unit = "cm", res = 1200)
  par(mfrow = c(1,1), mar = c(3, 3, 2, 2), mgp = c(1.7, 0.7, 0), bty = "n")
  plot(res_pr,'trace')
  dev.off()
  
  # ggplot(gene_dist_df, #%>% dplyr::filter(Gene == "group_2797"),
  #        aes(x = record,
  #            colour = Presence,
  #            fill = Presence)) +
  #   geom_bar() +
  #   facet_wrap(~SC) +
  #   theme_bw()
  
  # MCMC analysis
  mcmc_result <- BactDating::as.mcmc.resBactDating(res_pr)
  
  # Calculating ESS & Acceptance Rate
  calc_ess <- ess_calculation(mcmc_result)
  write.csv(calc_ess, "outputs/genomics/choosen_n703/calc_ess.csv", row.names = TRUE)
  
  # Figures! (still failed, margin error)
  png("pictures/genomics/choosen_n703/tree_trace2.png", width = 24, height = 12, unit = "cm", res = 1200)
  par(mfrow = c(1,1), mar = c(3, 3, 2, 2), mgp = c(1.7, 0.7, 0), bty = "n")
  pmcmc_trace(mcmc_result)
  dev.off()
  
  Sys.sleep(10) # wait 10 secs
  
  # Rooted tree!
  rooted_tree <- BactDating::initRoot(tre,d[,1]) # Incompatible dimensions because of d as matrix of (74,2)
  saveRDS(rooted_tree, "outputs/genomics/choosen_n703/rooted_tree.rds")
  
  png("pictures/genomics/choosen_n703/tree_rootedtree.png", width = 24, height = 12, unit = "cm", res = 1200)
  par(mfrow = c(1,1), mar = c(3, 3, 2, 2), mgp = c(1.7, 0.7, 0), bty = "n")
  plot(rooted_tree, show.tip.label = F)
  dev.off()
  
  Sys.sleep(10) # wait 10 secs
  
  png("pictures/genomics/choosen_n703/tree_roottotip.png", width = 24, height = 12, unit = "cm", res = 1200)
  par(mfrow = c(1,1), mar = c(3, 3, 2, 2), mgp = c(1.7, 0.7, 0), bty = "n")
  res_roottotip <- BactDating::roottotip(rooted_tree,d[,1])
  saveRDS(res_roottotip, "outputs/genomics/choosen_n703/res_roottotip.rds")
  dev.off()
  
}




# Some info about model selection:
# https://taming-the-beast.org/tutorials/NS-tutorial/

# Further analysis
# https://xavierdidelot.github.io/BactDating/articles/yourData.html
# https://xavierdidelot.github.io/BactDating/articles/Staph.html

# 
# modell <- c("mixedgamma", "strictgamma")
# res_post <- list()
# 
# for (i in 1:length(modell)){
#   
#   res_post[[i]]=bactdate(rooted,d,nbIts=1e6,
#                          initMu = effectiveSize(mcmc)["mu"],
#                          initAlpha = effectiveSize(mcmc)["alpha"],
#                          initSigma = effectiveSize(mcmc)["sigma"],
#                          model = modell[i], # try "strictgamma"
#                          showProgress = T)
#   
#   for (j in 1:length(res_post[[i]])){
#     
#     # MCMC traces
#     png(file = paste(modell[i], "_trace.png", sep = ""),
#         width = 1200, height = 650)
#     plot(res_post[[i]], 'trace')
#     dev.off()
#     
#     # treeCI
#     png(file = paste(modell[i], "_treeCI.png", sep = ""),
#         width = 1370, height = 960)
#     plot(res_post[[i]],'treeCI',show.tip.label = T)
#     dev.off()
#     
#     # treeRoot
#     png(file = paste(modell[i], "_treeRoot.png", sep = ""),
#         width = 1200, height = 910)
#     plot(res_post[[i]],'treeRoot',show.tip.label=T)
#     dev.off()
#     
#     # axisPhylo
#     png(file = paste(modell[i], "_axisPhylo.png", sep = ""),
#         width = 1300, height = 700)
#     out=extractSample(res_post[[i]],6)
#     par(mfrow=c(2,3),mar=c(5,5,0.5,0.5))
#     for (k in 1:6) {
#       plot(out[[k]])
#       axisPhylo(1,backward = F)
#     }
#     dev.off()
#   }
# }
# 
# 
# par(mfrow=c(1,1))
# tree=simcoaltree(2001:2015)
# plot(tree)
# ape::axisPhylo(backward=F)
# 
# obsphy=simobsphy(tree)
# r=roottotip(obsphy,2001:2015)
# 
# res2=bactdate(obsphy,2001:2015)
# plot(res2,'treeCI')