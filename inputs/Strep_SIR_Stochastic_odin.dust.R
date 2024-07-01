
library(odin.dust)
gen_sir <- odin.dust::odin_dust("inputs/sir_stochastic.R")

# Running the SIR model with dust
pars <- list(log_A_ini = (-3.77931571203353), # S_ini*10^(-5.69897) = 120 people; change A_ini into log10(A_ini)
             time_shift = 0.352271704195464,
             beta_0 = 0.0637307345664443,
             beta_1 = 0.174886960992199,
             beta_2 = 0.2,
             scaled_wane = (0.098), # 5.81837298310795E-05 for SIR model
             log_delta = (-4.55125666926139), # will be fitted to logN(-10, 0.7)
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
beta <- pars$beta_0 + pars$beta_0*pars$beta_1*sin(2*pi*((time_shift)+time)/365) + pars$beta_0*pars$beta_2*sin(2*pi*((time_shift*365)+time)/365)
# beta <- pars$beta_0*(1+pars$beta_1*sin(2*pi*(time_shift+time)/365))
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
## 1.1. Daily incidence
incidence <- read.csv("inputs/incidence.csv")

par(mfrow = c(1,1), mar = c(5.1, 5.1, 0.5, 0.5), mgp = c(3.5, 1, 0), las = 1)
cols <- c(S = "#8c8cd9", A = "darkred", D = "orange", R = "#999966", n_AD_daily = "#cc0099", n_AD_cumul = "green")
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

## 1.2. Weekly incidence
incidence_weekly <- read.csv("inputs/incidence_weekly.csv")
data_weekly <- incidence_weekly %>% 
  dplyr::group_by(weeks) %>% 
  dplyr::summarise(cases_weekly = sum(cases))

# Data preparation for model
min_date <- "2003-01-01"
max_date <- "2015-12-28"
all_date <- data.frame(allDate = seq.Date(from = as.Date(min_date),
                                          to = as.Date(max_date), 
                                          by = 1))
daily_incidence_modelled <- as.data.frame(t(x[5, , ]))
model_binds <- dplyr::bind_cols(all_date, daily_incidence_modelled)
model_binds$weeks <- cut(model_binds[,"allDate"], breaks="week")
model_weekly <- model_binds %>% 
  dplyr::group_by(weeks) %>% 
  dplyr::summarise(model_weekly = sum(V1))

data_plus_model <- dplyr::full_join(data_weekly, model_weekly,
                             by = c("weeks"))

# Plot!
png("pictures/data_plus_model.png", width = 17, height = 12, unit = "cm", res = 1200)
col_imD_weekly <- c(cases_weekly = "deepskyblue3",
                    model_weekly = "maroon")
ggplot(data_plus_model, aes(as.Date(weeks))) +
  geom_line(aes(y = cases_weekly, colour = "cases_weekly")) +
  geom_line(aes(y = model_weekly, colour = "model_weekly")) +
  scale_x_date() +
  scale_color_manual(values = col_imD_weekly,
                     name = "Cases",
                     breaks = c("cases_weekly", "model_weekly"),
                     labels = c("Data", "Model")
  ) +
  ggtitle("The Comparison of Model Output and Counts of Serotype 1 in England") +
  xlab("Year") +
  ylab("Serotype 1 Cases (Aggregated by Week)")
dev.off()

# Toy data creation ############################################################
# glimpse(x)
particle_1Data <- as.data.frame(x[5,,])
transposed_particle_1Data <- t(particle_1Data)

incidence_1particle <- tibble(day = 1:4745) %>% 
  bind_cols(transposed_particle_1Data)
# write.csv(incidence_1particle, file="outputs/SIS_daily_incidence_15particles.csv", row.names =F)
# write.csv(incidence_1particle, file="outputs/SIR_daily_incidence_15particles.csv", row.names =F)





# all_Data <- as.data.frame(x[5, 1:15, 1:4745])
# nrows <- length(all_Data[,,1]) # nrows = 4745
# ncols <- length(all_Data[,1,]) # ncols = 15
# 
# modified <- as.data.frame(matrix(all_Data, nrow = 4745, ncol = 15, byrow = T))
# glimpse(modified)
# 
# incidence_particles <- tibble(day = 1:4745) %>% 
#   bind_cols(modified)
# 
# write.csv(incidence_particles, file="outputs/SIS_daily_incidence_15particles.csv", row.names =F)
