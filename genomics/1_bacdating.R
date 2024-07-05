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
  
  tre <- BactDating::loadGubbins("raw_data/gubbins/n739_")
  
  # Temporary analysis coz' they don't share the dates (yet) ###################
  # tre <- loadGubbins("raw_data/gubbins/n739_")
  # datess <- as.vector(seq.Date(from = as.Date("2003-01-01"),
  #                              to = as.Date("2015-12-28"),
  #                              by = "year"))
  # truetreee=simcoaltree(datess)
  # tre=unroot(simobsphy(truetreee))
  # tre$unrec=runif(length(tre$edge.length))
  # plot(tre)
  # axisPhylo(backward = F)
  # edgelabels(round(100*tre$unrec))
  
  
  ##############################################################################
  tre_names <- as.data.frame(tre$tip.label) #%>% 
  # rename(ID = 'tre$tip.label')
  tre_names$ID <- substr(tre$tip.label, 1, 8)
  tre_names <- dplyr::left_join(tre_names, data, by = c("ID" = "ngsid"))
  tre_names <- dplyr::left_join(tre_names, dat, by = c("ID.y" = "ID"))
  tre_names <- tre_names %>% 
    dplyr::mutate(Earliestspecimendate = as.Date(Earliestspecimendate))
  
  # 2. BacDating ###############################################################
  # https://xavierdidelot.github.io/BactDating/articles/yourData.html
  d <- cbind(tre_names$Earliestspecimendate, tre_names$Earliestspecimendate+1) # 2-d matrix according to the articles above
  
  
  set.seed(0)
  
  res_pr <- BactDating::bactdate(tre,d,nbIts=nbIts, # Put 1e6 or 1e10 on hipercow
                                 model = "arc",
                                 showProgress = T)
  
  dir.create("outputs/genomics", FALSE, TRUE)
  saveRDS(res_pr, "outputs/genomics/mcmc_bacdating.rds")
  
  # Figures!
  fig <- plot(res_pr,'treeCI',show.tip.label = F)
  fig <- plot(res_pr,'trace')
  
  # MCMC analysis
  mcmc_result <- BactDating::as.mcmc.resBactDating(res_pr)
  
  # Calculating ESS & Acceptance Rate
  calc_ess <- ess_calculation(mcmc_result)
  write.csv(calc_ess, "outputs/genomics/calc_ess.csv", row.names = TRUE)
  
  # Figures! (still failed, margin error)
  fig <- pmcmc_trace(mcmc_result)
  
  Sys.sleep(10) # wait 10 secs
  
  # Rooted tree!
  rooted_tree <- BactDating::initRoot(tre,d[,1]) # Incompatible dimensions because of d as matrix of (74,2)
  saveRDS(rooted_tree, "outputs/genomics/rooted_tree.rds")
  
  Sys.sleep(10) # wait 10 secs
  
  res_roottotip <- BactDating::roottotip(rooted_tree,d[,1])
  saveRDS(res_roottotip, "outputs/genomics/res_roottotip.rds")
  
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