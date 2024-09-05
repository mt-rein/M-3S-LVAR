#### This script performs the simulation and saves the output ####

#### load packages and functions ####
library(conflicted)
library(dplyr)
library(flock)
library(lavaan)
library(MASS)
library(mcclust)
library(OpenMx)
library(purrr)
library(tidyr)
library(truncnorm)
library(RcppAlgos)
library(stringr)



# packages for parallel processing:
library(parabar)
library(parallel)

# load functions
source("step1.R")
source("step2.R")
source("step3.R")
source("do_sim.R")
source("auxiliary_functions.R")

#### define condition grid
cond <- expand.grid(replication = 1:50,
                    n = c(48, 96),
                    obs = c(25, 50, 100),
                    n_k = c(2, 4),
                    k_size = c("balanced", "unbalanced"),
                    rho_gen = c("low", "high"),
                    cluster_separation = c("low", "high"),
                    innovars = c("invariant", "random"))

# add seeds:
set.seed(123)
cond$seed <- sample(1:nrow(cond)*5, size = nrow(cond), replace = FALSE)
# add iteration number:
cond$iteration <- 1:nrow(cond)                                                  # (unique) number of each iteration


#### set up parallel computing ####
## open cluster
numCores <- parallel::detectCores() - 1
backend <- start_backend(numCores, cluster_type = "psock", backend_type = "async")


## load libraries in cluster
parabar::evaluate(backend, {
  library(conflicted)
  library(dplyr)
  library(flock)
  library(lavaan)
  library(MASS)
  library(mcclust)
  library(OpenMx)
  library(purrr)
  library(tidyr)
  library(truncnorm)
  library(RcppAlgos)
  library(stringr)
  
  source("step1.R")
  source("step2.R")
  source("step3_OpenMx_v1.R")
  source("do_sim.R")
  source("auxiliary_functions.R")
})

## load objects in cluster
export(backend, "cond")

#### perform simulation ####
start  <- Sys.time()
output <- par_lapply(backend, 1:nrow(cond), do_sim, cond = cond, outputfile = "test.csv", verbose = FALSE)
end <- Sys.time()

# close cluster:
stop_backend(backend)