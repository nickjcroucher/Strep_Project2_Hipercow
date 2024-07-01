# 2. Data Fitting ##############################################################
library(mcstate)
library(coda)
library(odin.dust)
library(dust)
library(GGally)
# library(socialmixr)

source("global/all_function.R") # Collected functions stored here!

# The anatomy of an mcstate particle filter, as noted above, consists of three main components: \n 
# 1. A set of observations to fit the model to, generated using mcstate::particle_filter_data(). \n 
# 2. A model to fit, which must be a dust generator, either dust::dust() or odin.dust::odin_dust(). \n 
# 3. A comparison function, which is an R function which calculates the likelihood of the state given the data at one time point.

# There is a calibration function in mcstate to fit our model to data.
# https://mrc-ide.github.io/mcstate/articles/sir_models.html

# To make my life easier I compile the Serotype 1 cases into a new object called sir_data
# data is fed as an input to mcstate::particle_filter_data
incidence <- read.csv("inputs/incidence.csv")

dt <- 1 # rate must be an integer; 0.25 to make it 4 days, I make it 1
sir_data <- mcstate::particle_filter_data(data = incidence,
                                          time = "day",
                                          rate = 1 / dt,
                                          initial_time = 0) # Initial time makes t0 start from 0 (not 1)

# Annotate the data so that it is suitable for the particle filter to use
rmarkdown::paged_table(sir_data)


## 2a. Model Load ##############################################################
# The model below is stochastic, closed system SADR model that I have created before
# I updated the code, filled the parameters with numbers;
# e.g.dt <- user(0) because if dt <- user() generates error during MCMC run
gen_sir <- odin.dust::odin_dust("inputs/sir_stochastic.R")

# This is part of sir odin model:
pars <- list(log_A_ini = (-5.69897), # S_ini*10^(log10(-5.69897)) = 120 people; change A_ini into log10(A_ini)
             time_shift = 0.2,
             beta_0 = 0.06565,
             beta_1 = 0.07, # in toy data the real value of beta_1 = 0.07
             max_wane = (-0.5),
             min_wane = (-4),
             scaled_wane = (0.5),
             log_delta = (-4.98),
             sigma_2 = 1
) # Serotype 1 is categorised to have the lowest carriage duration

# https://mrc-ide.github.io/odin-dust-tutorial/mcstate.html#/the-model-over-time
n_particles <- 50 # Trial n_particles = 50
filter <- mcstate::particle_filter$new(data = sir_data,
                                       model = gen_sir, # Use odin.dust input
                                       n_particles = n_particles,
                                       compare = case_compare,
                                       seed = 1L)

filter$run(pars)

# Variance and particles estimation (as suggested by Rich)
# parallel::detectCores() # I have 4 cores
# x <- replicate(30, filter$run(pars))
# x
# 
# var(x)
# # [1] 266.5598
# 69 / 267 # Trial 69 particles for 267 var; how many particles are needed?
# # [1] 0.258427
# 69  / 4 # change by the factor of 4
# # [1] 17.25
# 69  / 4  /4
# # [1] 4.3125
# 69  / 4  /4 /4
# # [1] 1.078125 # so the factor is 4*4*4 to finally get roughly 1 particle
# 4 * 4 * 4
# # [1] 64
# 4 * 4 * 4 * 500
# # [1] 32000 --> particles needed for var(x) = 1

# Update n_particles based on calculation in 4 cores with var(x) ~ 267: 32000

priors <- prepare_priors(pars)
proposal_matrix <- diag(300, 6) # previously 200
proposal_matrix <- (proposal_matrix + t(proposal_matrix)) / 2
rownames(proposal_matrix) <- c("log_A_ini", "time_shift", "beta_0", "beta_1", "scaled_wane", "log_delta")
colnames(proposal_matrix) <- c("log_A_ini", "time_shift", "beta_0", "beta_1", "scaled_wane", "log_delta")

mcmc_pars <- prepare_parameters(initial_pars = pars, priors = priors, proposal = proposal_matrix, transform = transform)

# n_steps <- 100 #1e6

# I change pmcmc_run into a function that involve control inside:
# pmcmc_run <- mcstate::pmcmc(mcmc_pars, filter_deterministic, control = control)

# Directory for saving the outputs
dir.create("outputs", FALSE, TRUE)

# Trial combine pMCMC + tuning #################################################
pmcmc_run_plus_tuning <- function(n_particles, n_steps){
  filter <- mcstate::particle_filter$new(data = sir_data,
                                         model = gen_sir, # Use odin.dust input
                                         n_particles = n_particles,
                                         compare = case_compare,
                                         seed = 1L)
  
  # Use deterministic model by add filter_deterministic
  # https://mrc-ide.github.io/mcstate/articles/deterministic.html
  # Index function is optional when only a small number of states are used in comparison function.
  filter_deterministic <- mcstate::particle_deterministic$new(data = sir_data,
                                                              model = gen_sir,
                                                              compare = case_compare
                                                              # index = index
  )
  
  
  control <- mcstate::pmcmc_control(n_steps = n_steps,
                                    rerun_every = 50,
                                    rerun_random = TRUE,
                                    progress = TRUE)
  
  # The pmcmc
  pmcmc_result <- mcstate::pmcmc(mcmc_pars, filter_deterministic, control = control)
  pmcmc_result
  saveRDS(pmcmc_result, "outputs/pmcmc_result.rds")
  
  new_proposal_mtx <- cov(pmcmc_result$pars)
  write.csv(new_proposal_mtx, "outputs/new_proposal_mtx.csv", row.names = TRUE)
  
  lpost_max <- which.max(pmcmc_result$probabilities[, "log_posterior"])
  write.csv(as.list(pmcmc_result$pars[lpost_max, ]),
            "outputs/initial.csv", row.names = FALSE)
  
  # Further processing for thinning chains
  mcmc1 <- pmcmc_further_process(n_steps, pmcmc_result)
  write.csv(mcmc1, "outputs/mcmc1.csv", row.names = TRUE)
  
  # Calculating ESS & Acceptance Rate
  calc_ess <- ess_calculation(mcmc1)
  write.csv(calc_ess, "outputs/calc_ess.csv", row.names = TRUE)
  
  # Figures! (still failed, margin error)
  fig <- pmcmc_trace(mcmc1)
  
  Sys.sleep(10) # wait 10 secs before conducting tuning
  
  # New proposal matrix
  new_proposal_matrix <- as.matrix(read.csv("outputs/new_proposal_mtx.csv"))
  new_proposal_matrix <- new_proposal_matrix[, -1]
  new_proposal_matrix <- apply(new_proposal_matrix, 2, as.numeric)
  new_proposal_matrix <- new_proposal_matrix/100 # Lilith's suggestion
  new_proposal_matrix <- (new_proposal_matrix + t(new_proposal_matrix)) / 2
  rownames(new_proposal_matrix) <- c("log_A_ini", "time_shift", "beta_0", "beta_1", "scaled_wane", "log_delta")
  colnames(new_proposal_matrix) <- c("log_A_ini", "time_shift", "beta_0", "beta_1", "scaled_wane", "log_delta")
  # isSymmetric(new_proposal_matrix)
  
  tune_mcmc_pars <- prepare_parameters(initial_pars = pars, priors = priors, proposal = new_proposal_matrix, transform = transform)
  
  # Including adaptive proposal control
  # https://mrc-ide.github.io/mcstate/reference/adaptive_proposal_control.html
  tune_control <- mcstate::pmcmc_control(n_steps = n_steps,
                                         n_chains = 4,
                                         rerun_every = 50,
                                         rerun_random = TRUE,
                                         progress = TRUE,
                                         adaptive_proposal = adaptive_proposal_control(initial_vcv_weight = 1000,
                                                                                       # initial_scaling = 1,
                                                                                       scaling_increment = NULL,
                                                                                       # log_scaling_update = T,
                                                                                       acceptance_target = 0.234,
                                                                                       forget_rate = 0.2,
                                                                                       forget_end = Inf,
                                                                                       adapt_end = Inf
                                         ))
  
  filter <- mcstate::particle_filter$new(data = sir_data,
                                         model = gen_sir, # Use odin.dust input
                                         n_particles = n_particles,
                                         compare = case_compare,
                                         seed = 1L
                                         )
  
  # The pmcmc
  tune_pmcmc_result <- mcstate::pmcmc(tune_mcmc_pars, filter_deterministic, control = tune_control)
  tune_pmcmc_result
  saveRDS(tune_pmcmc_result, "outputs/tune_pmcmc_result.rds")
  
  tune_lpost_max <- which.max(tune_pmcmc_result$probabilities[, "log_posterior"])
  mcmc_lo_CI <- apply(tune_pmcmc_result$pars, 2, function(x) quantile(x, probs = 0.025))
  mcmc_hi_CI <- apply(tune_pmcmc_result$pars, 2, function(x) quantile(x, probs = 0.975))
  
  binds_tune_initial <- rbind(as.list(tune_pmcmc_result$pars[tune_lpost_max, ]), mcmc_lo_CI, mcmc_hi_CI)
  t_tune_initial <- t(binds_tune_initial)
  colnames(t_tune_initial) <- c("values", "low_CI", "high_CI")
  
  write.csv(t_tune_initial,
            "outputs/tune_initial_with_CI.csv", row.names = FALSE)
  
  # Further processing for thinning chains
  mcmc2 <- tuning_pmcmc_further_process(n_steps, tune_pmcmc_result)
  # mcmc2 <- coda::as.mcmc(cbind(
  #   tune_pmcmc_result$probabilities, tune_pmcmc_result$pars))
  write.csv(mcmc2, "outputs/mcmc2.csv", row.names = TRUE)
  
  # Calculating ESS & Acceptance Rate
  tune_calc_ess <- ess_calculation(mcmc2)
  write.csv(tune_calc_ess, "outputs/tune_calc_ess.csv", row.names = TRUE)
  
  # Figures! (still failed, margin error)
  fig <- pmcmc_trace(mcmc2)
  
  ##############################################################################
  # MCMC Diagnostics
  
  # 1. Gelman-Rubin Diagnostic
  # https://cran.r-project.org/web/packages/coda/coda.pdf
  # png("pictures/diag_gelman_rubin.png", width = 17, height = 12, unit = "cm", res = 1200)
  figs_gelman_init <- diag_init_gelman_rubin(tune_pmcmc_result)
  fig <- diag_cov_mtx(figs_gelman_init)
  fig <- diag_gelman_rubin(figs_gelman_init)
  # dev.off()
  
  # 2. Autocorrelation
  # png("pictures/diag_aucorr.png", width = 17, height = 12, unit = "cm", res = 1200)
  fig <- diag_aucorr(mcmc2)
  # dev.off()
  
  # png("pictures/diag_ggpairs.png", width = 17, height = 12, unit = "cm", res = 1200)
  fig <- GGally::ggpairs(as.data.frame(tune_pmcmc_result$pars))
  # dev.off()
  
}

# pmcmc_run_plus_tuning(40000, 1e3)
