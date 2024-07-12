# Libary
library(odin.dust)

# Set up simulation parameters
n_particles <- 10L
n_times <- 1000L

# Compile model
gen_sir <- odin.dust::odin_dust("./inputs/sir_stochastic.R")

# Run model
sir_model <- gen_sir$new(pars = list(log_A_ini = (-5.69897), # S_ini*10^(log10(-5.69897)) = 120 people; change A_ini into log10(A_ini)
                                     time_shift_1 = 0.2,
                                     time_shift_2 = 0.2,
                                     beta_0 = 0.06565,
                                     beta_1 = 0.07, # in toy data the real value of beta_1 = 0.07
                                     beta_2 = 0.2,
                                     max_wane = (-0.5),
                                     min_wane = (-5),
                                     scaled_wane = (0.01),
                                     log_delta = (-4.98),
                                     sigma_2 = 1
                          ),
                         time = 1L,
                         n_particles = n_particles,
                         n_threads = 4L,
                         seed = 1L)


# Plot output
x <- array(NA, dim = c(sir_model$info()$len, n_particles, n_times))

for (t in seq_len(n_times)) {
  x[ , , t] <- sir_model$run(t)
}
time <- x[1, 1, ]
x <- x[-1, , ]

par(mar = c(4.1, 5.1, 0.5, 0.5), las = 1)
cols <- c(S = "#8c8cd9", A = "#cc0044", D = "#999966", R = "orange")
matplot(time, t(x[1, , ]), type = "l",
        xlab = "Time", ylab = "Number of individuals",
        col = cols[["S"]], lty = 1, ylim = range(x))
matlines(time, t(x[2, , ]), col = cols[["A"]], lty = 1)
matlines(time, t(x[3, , ]), col = cols[["D"]], lty = 1)
matlines(time, t(x[4, , ]), col = cols[["R"]], lty = 1)
legend("left", lwd = 1, col = cols, legend = names(cols), bty = "n")
