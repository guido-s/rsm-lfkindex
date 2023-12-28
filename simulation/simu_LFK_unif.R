## R packages meta and metasens called internally
##
library("parallel")
library("tictoc")
library("dplyr")


## Loading functions
##
source("funcs_simu.R")

## set seed
set.seed(1919)


## Set cores
##
numCores <- detectCores()
numCores <- max(numCores - 2, 1)
numCores


## 10000 iterations
##
nsim <- 10000


## Sample size distribution
##
n <- draw_n_unif(nsim)
n <- n$e + n$c


## Define scenarios (rows)
##
param.grid <-
  expand.grid( # all combinations
    N_studies = c(10, 20, 50, 100),
    theta_treat = c(0, 1, 2),
    sd_treat = c(1, 1.5, 2),
    theta_control = 2,
    sd_control = 1,
    tau = c(0, 0.25, 0.5),
    rho = c(0, -0.5, -0.9),
    p_small = 0.1,
    p_large = 0.9,
    n_small = round(quantile(n, 0.1)),
    n_large = round(quantile(n, 0.9))
  )
##
attr(param.grid, "out.attrs") <- NULL
param.grid <- cbind(scenario = seq_len(nrow(param.grid)), param.grid)
##
param.names <- names(param.grid)


## Run simulation
##
params <- list()
##
for (i in seq_len(nrow(param.grid)))
  params[[i]] <- as.list(param.grid[i, ])
##
tic()
sim.res <-
  mclapply(params, runsim, nsim = nsim, mc.cores = numCores,
           func.n = draw_n_unif)
toc()

save(sim.res, param.grid, param.names, nsim, params,
     file = "simu_LFK_unif.rda")
