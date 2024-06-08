
case_compare <- function(state, observed, pars = NULL) {
  exp_noise <- 1e4
  
  incidence_modelled <- state[6, , drop = TRUE] # (incidence based on model's "n_AD_daily" from gen_sir)
  incidence_observed <- observed$cases # daily new cases
  lamb <- incidence_modelled +
    rexp(n = length(incidence_modelled), rate = exp_noise)
  dpois(x = incidence_observed, lambda = lamb, log = TRUE)
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
    list(mcstate::pmcmc_parameter("time_shift", 0.2, min = 0, max = 1,
                                  prior = function(s) dunif(s, min = 0, max = 1, log = TRUE)), # ~Uniform[0,1] in the proportion of 365 days
         mcstate::pmcmc_parameter("beta_0", 0.06565, min = 0, max = 0.8,
                                  prior = function(s) dgamma(s, shape = 1, scale = 0.1, log = TRUE)), # draws from gamma distribution dgamma(1, 0.2) --> exp dist)
         mcstate::pmcmc_parameter("beta_1", 0.07, min = 0, max = 0.8,
                                  prior = function(s) dgamma(s, shape = 1, scale = 0.1, log = TRUE)), # draws from gamma distribution dgamma(1, 0.2) --> exp dist
         mcstate::pmcmc_parameter("wane", 0.002, min = 0, max = 0.8,
                                  prior = function(s) dgamma(s, shape = 1, scale = 0.1, log = TRUE)), # draws from gamma distribution dgamma(1, 0.2) --> exp dist
         mcstate::pmcmc_parameter("log_delta", (-4.7), min = (-10), max = 0.7,
                                  prior = function(s) dunif(s, min = (-10), max = 0.7, log = TRUE)), # logN distribution for children & adults (Lochen et al., 2022)
         mcstate::pmcmc_parameter("sigma_2", 1, min = 0, max = 10,
                                  prior = function(s) dgamma(s, shape = 1, scale = 1, log = TRUE)) # shape = 1 , scale = 1 to capture 5 days/more (or dgamma(2.5, 0.5))?
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
