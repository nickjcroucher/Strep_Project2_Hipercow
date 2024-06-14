
case_compare <- function(state, observed, pars = NULL) {
  exp_noise <- 1e4
  
  # incidence based on model's "n_SI_daily" from gen_sir$new(pars = list(), time = 0, n_particles = 1L)$info()
  incidence_modelled <- state[5, , drop = TRUE] # n_SI_daily is located at state[5, , ]
  
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
  time_shift <- pars[["time_shift"]]
  beta_0 <- pars[["beta_0"]]
  beta_1 <- pars[["beta_1"]]
  wane <- pars[["wane"]]
  log_delta <- pars[["log_delta"]]
  sigma_2 <- pars[["sigma_2"]]
  
  list(time_shift = time_shift,
       beta_0 = beta_0,
       beta_1 = beta_1,
       wane = wane,
       log_delta = log_delta,
       sigma_2 = sigma_2)
  
}

transform <- function(pars) {
  parameter_transform(pars)
}

prepare_parameters <- function(initial_pars, priors, proposal, transform) {
  
  mcmc_pars <- mcstate::pmcmc_parameters$new(
    list(mcstate::pmcmc_parameter("time_shift", 0.1, min = 0, max = 1,
                                  prior = priors$time_shift), # ~Uniform[0,1] in the proportion of 365 days
         mcstate::pmcmc_parameter("beta_0", 0.16565, min = 0, max = 0.8,
                                  prior = priors$beta_0), # draws from gamma distribution dgamma(1, 0.2) --> exp dist)
         mcstate::pmcmc_parameter("beta_1", 0.05, min = 0, max = 0.8,
                                  prior = priors$beta_1), # draws from gamma distribution dgamma(1, 0.2) --> exp dist
         mcstate::pmcmc_parameter("wane", 0.002, min = 0, max = 0.8,
                                  prior = priors$wane), # draws from gamma distribution dgamma(1, 0.2) --> exp dist
         mcstate::pmcmc_parameter("log_delta", (-4.7), min = (-10), max = 0.7,
                                  prior = priors$log_delta), # logN distribution for children & adults (Lochen et al., 2022)
         mcstate::pmcmc_parameter("sigma_2", 1, min = 0, max = 10,
                                  prior = priors$sigma_2) # shape = 1 , scale = 1 to capture 5 days/more (or dgamma(2.5, 0.5))?
    ),
    proposal = proposal,
    transform = transform)
  
}

prepare_priors <- function(pars) {
  priors <- list()
  
  priors$time_shift <- function(s) {
    dunif(s, min = 0, max = 1, log = TRUE)
  }
  priors$beta_0 <- function(s) {
    dgamma(s, shape = 1, scale = 0.1, log = TRUE)
  }
  priors$beta_1 <- function(s) {
    dgamma(s, shape = 1, scale = 0.1, log = TRUE)
  }
  priors$wane <- function(s) {
    dgamma(s, shape = 1, scale = 0.1, log = TRUE)
  }
  priors$log_delta <- function(s) {
    dunif(s, min = (-10), max = 0.7, log = TRUE)
  }
  priors$sigma_2 <- function(s) {
    dgamma(s, shape = 1, scale = 1, log = TRUE)
  }
  priors
}

pmcmc_further_process <- function(n_steps, pmcmc_result) {
  processed_chains <- mcstate::pmcmc_thin(pmcmc_result, burnin = n_steps/2, thin = 2)
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
  processed_chains <- mcstate::pmcmc_thin(tune_pmcmc_result, burnin = n_steps/2, thin = 2)
  parameter_mean_hpd <- apply(processed_chains$pars, 2, mean)
  parameter_mean_hpd
  
  mcmc_tuning_result <- coda::as.mcmc(cbind(processed_chains$probabilities, processed_chains$pars))
  mcmc_tuning_result
}
