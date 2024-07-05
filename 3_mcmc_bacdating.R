
# The pmcmc run by using Hipercow! #############################################
# https://mrc-ide.github.io/hipercow/articles/windows.html
library(hipercow)
library(parallel)
library(ids)
# sudo mount -a

# See filesystems and paths, CHANGE wd to those in /etc/fstab
# DO NOT CHANGE THE ARRANGEMENT OF THESE COMMANDS!!!

hipercow_init(driver = "windows")
hipercow_configure("windows", r_version = "4.4.0")
windows_authenticate() # authenticate by using DIDE account
windows_check()
hipercow_configuration() # for troubleshooting
# hipercow_hello() # test job

hipercow_environment_create(sources = "genomics/1_bacdating.R")
hipercow_provision()

# Check the installed packages again by using hipercow_configuration()
# hipercow_configuration()

# If automatic install failed (don't know why), use pkgdepends.txt!
# install.packages("pkgdepends")
# hipercow_provision()

# https://mrc-ide.github.io/hipercow/reference/hipercow_resources.html
resources <- hipercow::hipercow_resources(cores = 20,
                                          max_runtime = "3d",
                                          memory_per_node = "64G",
)


# Now pmcmc_run is a function:
# pmcmc_run <- function(n_particles, n_steps)
run_tree <- task_create_expr(run_bacdating(1e6), # try 1e8
                             resources = resources
)

# Something related to test the submitted job
task_status(run_tree)
task_result(run_tree)
task_log_show(run_tree)
task_info(run_tree)
task_info(run_tree)$times

