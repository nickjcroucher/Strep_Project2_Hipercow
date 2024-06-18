
library(odin.dust)
gen_sir <- odin.dust::odin_dust("inputs/sir_stochastic.R")

# Running the SIR model with dust
pars <- list(log_A_ini = (-5.69897), # S_ini*10^(log10(-5.69897)) = 120 people; change A_ini into log10(A_ini)
             time_shift = 0.2,
             beta_0 = 0.06565,
             beta_1 = 0.07,
             wane = 0.002,
             log_delta = (-4.98), # will be fitted to logN(-10, 0.7)
             sigma_1 = (1/15.75),
             sigma_2 = (1)
) # Serotype 1 is categorised to have the lowest carriage duration

sir_model <- gen_sir$new(pars = pars,
                         time = 1,
                         n_particles = 15L,
                         n_threads = 4L,
                         seed = 1L)

# sir_model$state() # test array OR matrix state

# update_state is required "every single time" to run & produce matrix output (don't know why)
sir_model$update_state(pars = pars,
                       time = 0) # make sure time is 0

# all_date <- incidence$day
# all_date <- data.frame(col = integer(4745))
n_times <- 4745 # 4745 or similar to the number of date range (of the provided data), or try 500 for trial
n_particles <- 15
x <- array(NA, dim = c(sir_model$info()$len, n_particles, n_times))

# Beta check
time <- seq(1, n_times, 1)
time_shift <- 70
beta <- pars$beta_0*(1+pars$beta_1*sin(2*pi*(time_shift+time)/365))
max(beta)
min(beta)

# R0 estimation (R0 changes due to seasonality)
R0 <- (beta/(pars$log_delta+pars$sigma_1)) +  ((pars$log_delta)*(beta)) / ((pars$log_delta + 192/(4064*4745))*(pars$sigma_2 + 192/(4064*4745))) # print R0
max(R0)
min(R0)
# plot(time, R0)
# pars$beta_1/(pars$delta) + (pars$qu*pars$delta)/(pars$delta*pars$sigma) # print R0


for (t in seq_len(n_times)) {
  x[ , , t] <- sir_model$run(t)
}
time <- x[1, 1, ] # because in the position of [1, 1, ] is time
x <- x[-1, , ] # compile all matrix into 1 huge df, delete time (position [-1, , ])
library(tidyverse)
glimpse(x)

## 1. Data Load ################################################################
incidence <- read.csv("inputs/incidence.csv")

par(mar = c(5.1, 5.1, 0.5, 0.5), mgp = c(3.5, 1, 0), las = 1)
cols <- c(S = "#8c8cd9", A = "darkred", D = "orange", R = "#999966", n_AD_daily = "#cc0099", n_AD_cumul = "green")
# matplot(time, t(x[1, , ]), type = "l",
#         xlab = "Time", ylab = "Number of individuals",
#         col = cols[["S"]], lty = 1, ylim = range(x))
matplot(time, t(x[5, , ]), type = "l",
        xlab = "Time", ylab = "Number of individuals",
        col = cols[["n_AD_daily"]], lty = 1)#, ylim = max(x[2,,]))

matlines(incidence$day, incidence$cases, type = "l", col = "steelblue")

# matlines(time, t(x[2, , ]), col = cols[["A"]], lty = 1)
# matlines(time, t(x[3, , ]), col = cols[["D"]], lty = 1)
# matlines(time, t(x[4, , ]), col = cols[["R"]], lty = 1)
# matlines(time, t(x[5, , ]), col = cols[["n_AD_daily"]], lty = 1)
# matlines(time, t(x[6, , ]), col = cols[["n_AD_cumul"]], lty = 1)
legend("left", lwd = 1, col = cols, legend = names(cols), bty = "n")
max(x[5,,]) # Check max n_AD_daily
max(x[3,,]) # Check max D

# Toy data creation ############################################################
# glimpse(x)
new_toyData <- as.data.frame(x[5, 1, 1:4745])
colnames(new_toyData) <- "cases"
glimpse(new_toyData)

incidence <- tibble(day = 1:4745) %>% 
  bind_cols(new_toyData)

# write.csv(incidence, file="inputs/incidence_toyData.csv", row.names =F)
