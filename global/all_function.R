
case_compare <- function(state, observed, pars = NULL) {
  exp_noise <- 1e4
  
  # incidence based on model's "n_AD_daily" from gen_sir$new(pars = list(), time = 0, n_particles = 1L)$info()
  incidence_modelled <- state[6, , drop = TRUE] # n_AD_daily is located at state[6, , ]
  
  # incidence based on data
  incidence_observed <- observed$cases # daily new cases
  
  n <- ncol(state)
  lamb <- incidence_modelled + rexp(n, exp_noise)
  loglik_cases <- dpois(x = incidence_observed, lambda = lamb, log = T)
  
  return(loglik_cases)
}

# That transform function
# https://github.com/mrc-ide/mcstate/blob/da9f79e4b5dd421fd2e26b8b3d55c78735a29c27/tests/testthat/test-if2.R#L40
# https://github.com/mrc-ide/mcstate/issues/184
parameter_transform <- function(pars) {
  log_A_ini <- pars[["log_A_ini"]]
  time_shift_1 <- pars[["time_shift_1"]]
  time_shift_2 <- pars[["time_shift_2"]]
  beta_0 <- pars[["beta_0"]]
  beta_1 <- pars[["beta_1"]]
  beta_2 <- pars[["beta_2"]]
  scaled_wane <- pars[["scaled_wane"]]
  log_delta <- pars[["log_delta"]]
  # sigma_2 <- pars[["sigma_2"]]

  list(log_A_ini = log_A_ini,
       time_shift_1 = time_shift_1,
       time_shift_2 = time_shift_2,
       beta_0 = beta_0,
       beta_1 = beta_1,
       beta_2 = beta_2,
       scaled_wane = scaled_wane,
       log_delta = log_delta#,
       # sigma_2 = sigma_2
       )
  
}

transform <- function(pars) {
  parameter_transform(pars)
}

prepare_parameters <- function(initial_pars, priors, proposal, transform) {
  
  mcmc_pars <- mcstate::pmcmc_parameters$new(
    list(mcstate::pmcmc_parameter("log_A_ini", (-4.69897), min = (-10), max = 0,
                                  prior = priors$log_A_ini),
         mcstate::pmcmc_parameter("time_shift_1", 0.2, min = 0, max = 1,
                                  prior = priors$time_shifts),
         mcstate::pmcmc_parameter("time_shift_2", 0.2, min = 0, max = 1,
                                  prior = priors$time_shifts),
         mcstate::pmcmc_parameter("beta_0", 0.06565, min = 0, max = 0.8,
                                  prior = priors$betas),
         mcstate::pmcmc_parameter("beta_1", 0.07, min = 0, max = 1,
                                  prior = priors$betas),
         mcstate::pmcmc_parameter("beta_2", 0.2, min = 0, max = 1,
                                  prior = priors$betas),
         mcstate::pmcmc_parameter("scaled_wane", (0.5), min = (0), max = 1,
                                  prior = priors$scaled_wane),
         mcstate::pmcmc_parameter("log_delta", (-4.98), min = (-10), max = 0.7,
                                  prior = priors$log_delta)#,
         # mcstate::pmcmc_parameter("sigma_2", 1, min = 0, max = 10,
         #                          prior = priors$sigma_2)
    ),
    proposal = proposal,
    transform = transform)
  
}

prepare_priors <- function(pars) {
  priors <- list()
  
  priors$log_A_ini <- function(s) {
    dunif(s, min = (-10), max = 0, log = TRUE)
  }
  priors$time_shifts <- function(s) {
    dunif(s, min = 0, max = 1, log = TRUE)
  }
  priors$betas <- function(s) {
    dgamma(s, shape = 1, scale = 0.1, log = TRUE)
  }
  priors$scaled_wane <- function(s) {
    dbeta(s, shape1 = 1.25, shape2 = 1.25, log = TRUE)
  }
  priors$log_delta <- function(s) {
    dunif(s, min = (-10), max = 0.7, log = TRUE)
  }
  # priors$sigma_2 <- function(s) {
  #   dgamma(s, shape = 1, scale = 1, log = TRUE)
  # }
  priors
}

pmcmc_further_process <- function(n_steps, pmcmc_result) {
  processed_chains <- mcstate::pmcmc_thin(pmcmc_result, burnin = n_steps*0.95, thin = NULL)
  parameter_mean_hpd <- apply(processed_chains$pars, 2, mean)
  parameter_mean_hpd
  
  mcmc1 <- coda::as.mcmc(cbind(pmcmc_result$probabilities, pmcmc_result$pars))
  mcmc1
}

ess_calculation <- function(mcmc1){
  calc <- list(ess = coda::effectiveSize(mcmc1),
               acceptance_rate = 1 - coda::rejectionRate(mcmc1))
  calc
}

pmcmc_trace <- function(mcmc1) {
  plot(mcmc1) # to save the figures into pdf
  # png("pictures/mcmc1.png", res = 1200)
  # plot(mcmc1)
  # dev.off()
}

################################################################################
# Tuning functions
tuning_pmcmc_further_process <- function(n_steps, tune_pmcmc_result) {
  processed_chains <- mcstate::pmcmc_thin(tune_pmcmc_result, burnin = n_steps*0.95, thin = 2)
  parameter_mean_hpd <- apply(processed_chains$pars, 2, mean)
  parameter_mean_hpd
  
  tune_pmcmc_result <- coda::as.mcmc(cbind(processed_chains$probabilities, processed_chains$pars))
  tune_pmcmc_result
}

################################################################################
# MCMC Diagnostics
# 1. Gelman-Rubin Diagnostic
# https://cran.r-project.org/web/packages/coda/coda.pdf

diag_init_gelman_rubin <- function(tune_pmcmc_result){
  n_chains <- 4 # tune_control$n_chains
  n_samples <- nrow(tune_pmcmc_result$pars)/n_chains
  
  # Split the parameter samples and probabilities by chains
  chains <- lapply(1:n_chains, function(i) {
    start <- (i - 1) * n_samples + 1
    end <- i * n_samples
    list(
      pars = tune_pmcmc_result$pars[start:end, ],
      probabilities = tune_pmcmc_result$probabilities[start:end, ]
    )
  })
  
  # Convert chains to mcmc objects
  mcmc_chains <- lapply(chains, function(chain) {
    as.mcmc(cbind(chain$probabilities, chain$pars))
  })
  
  # Combine chains into a list
  mcmc_chains_list <- do.call(mcmc.list, mcmc_chains)
  mcmc_chains_list
}

diag_cov_mtx <- function(mcmc_chains_list) {
  # print("Covariance matrix of mcmc2")
  cov(as.matrix(mcmc_chains_list))
}

diag_gelman_rubin <- function(mcmc_chains_list) {
  # print("Gelman-Rubin diagnostic")
  gelman_plot <- coda::gelman.plot(mcmc_chains_list,
                    bin.width = 10,
                    max.bins = 50,
                    confidence = 0.95,
                    transform = FALSE,
                    autoburnin=TRUE,
                    auto.layout = TRUE)
  # ask, col, lty, xlab, ylab, type, ...)
  
  coda::gelman.diag(mcmc_chains_list,
                    confidence = 0.95,
                    transform=FALSE,
                    autoburnin=TRUE,
                    multivariate=F) # Change multivariate = F instead of T
  
}

# 2. Autocorrelation plots
diag_aucorr <- function(mcmc2){
  for (name in colnames(mcmc2)){
    print(coda::acfplot(mcmc2[, name], main = name))
  }
}



