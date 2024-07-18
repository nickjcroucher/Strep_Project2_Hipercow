# Libary
library(odin.dust)

# Set up simulation parameters
n_particles <- 10L
n_times <- 4745L

# Compile model
gen_sir <- odin.dust::odin_dust("./inputs/sir_stochastic.R")

# Run model
sir_model <- gen_sir$new(pars = list(log_A_ini = (-3.776898385), # S_ini*10^(log10(-5.69897)) = 120 people; change A_ini into log10(A_ini)
                                     time_shift_1 = 0.111895853,
                                     time_shift_2 = 0.344809421,
                                     beta_0 = 0.031878385,
                                     beta_1 = 0.147861631, # in toy data the real value of beta_1 = 0.07
                                     beta_2 = 0.197702508,
                                     max_wane = (-0.5),
                                     min_wane = (-5),
                                     scaled_wane = (0),
                                     log_delta = (-4.554225071),
                                     sigma_2 = 1
                          ),
                         time = 1L,
                         n_particles = n_particles,
                         n_threads = 4L,
                         seed = 1L)


# Plot output
simulation_df <-
  reshape2::melt(x) %>%
  dplyr::rename(
    Index = Var1,
    Replicate = Var2,
    Time = Var3
  ) %>%
  dplyr::filter(Index < 5) %>%
  dplyr::mutate(
    Compartment = 
      dplyr::case_when(
        Index == 1 ~ "Asymptomatic",
        Index == 2 ~ "Diseased",
        Index == 3 ~ "Susceptible",
        Index == 4 ~ "Recovered"
      )
  )

ggplot(simulation_df,
       aes(x = Time,
           y = value,
           colour = Compartment,
           group = interaction(Compartment,Replicate))) +
  geom_line() +
  scale_y_continuous(trans = "log1p") +
  theme_bw()

