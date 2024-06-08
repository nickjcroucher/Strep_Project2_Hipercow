
# The pmcmc run by using Hipercow! #############################################
# https://mrc-ide.github.io/hipercow/articles/windows.html
library(hipercow)
# sudo mount -a

# See filesystems and paths, CHANGE wd to those in /etc/fstab
# DO NOT CHANGE THE ARRANGEMENT OF THESE LINES BELOW!!!

hipercow_init(driver = "windows")
hipercow_configure("windows", r_version = "4.4.0")
windows_authenticate() # authenticate by using DIDE account
windows_check()
hipercow_configuration() # for troubleshooting
# hipercow_hello() # test job

hipercow_environment_create(name = "mcState_Model",
                            sources = c("1_data_preparation.R", "2_pmcmc.R"),
                            globals = "global/all_function.R",
                            overwrite = T,
                            check = T, # check if error occurs
                            root = "/home/ron/net/home/Strep_Project2_Hipercow"
                            )
hipercow_provision(environment = "mcState_Model",
                   root = "/home/ron/net/home/Strep_Project2_Hipercow")

# Check the installed packages again by using hipercow_configuration()
hipercow_configuration()

# If automatic install failed (don't know why), use pkgdepends.txt!
# install.packages("pkgdepends")
# hipercow_provision()

# https://mrc-ide.github.io/hipercow/reference/hipercow_resources.html
resources <- hipercow::hipercow_resources(cores = 20,
                                          max_runtime = "3d",
                                          memory_per_node = "64G",
                                          )

# https://mrc-ide.github.io/hipercow/reference/hipercow_parallel.html
parallel <- hipercow::hipercow_parallel(method = "parallel",
                                        cores_per_process = 20,
                                        environment = "mcState_Model")

# Now pmcmc_run is a function:
# pmcmc_run <- function(n_particles, n_steps)
id_single <- task_create_expr(pmcmc_run(40000, 1e6), # Update n_particles = 32000, n_steps = 1e6?
                       resources = resources
                       )

# Something related to test the submitted job
task_status(id_single)
task_result(id_single)
task_log_show(id_single)
task_info(id_single)

# Trial parallel job submission:
id_parallel <- task_create_expr(
  parallel::clusterApply(NULL, 1:20, pmcmc_run(500, 100),
    c(Sys.getpid(), hipercow_parallel_get_cores()),
    parallel = hipercow_parallel("parallel"),
    resources = resources))

# Something related to test the submitted job
task_status(id_parallel)
task_result(id_parallel)
task_log_show(id_parallel)
task_info(id_parallel)

# Something related to test the submitted job
# id <- task_create_expr(sessionInfo())
task_status(id)
task_result(id)
task_log_show(id)
task_info(id)
