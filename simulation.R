#### This script performs the simulation and saves the output ####

#### load packages and functions ####
library(conflicted)
library(dplyr)
library(flock)
library(lavaan)
library(MASS)
library(mcclust)
library(purrr)
library(tidyr)
library(truncnorm)



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
cond <- expand.grid(replication = 1,
                    n = c(48, 96),
                    obs = c(25, 50, 75, 100),
                    n_k = c(2, 4),
                    k_size = c("balanced", "unbalanced"),
                    rho_gen = c("small", "medium", "large", "very large"),
                    similarity = c("similar", "dissimilar"),
                    innovars = c("equal", "random"))

# add seeds:
set.seed(123)
cond$seed <- sample(1:nrow(cond)*5, size = nrow(cond), replace = FALSE)
# add iteration number:
cond$iteration <- 1:nrow(cond)                                                  # (unique) number of each iteration


# #### set up parallel computing ####
# ## open cluster
# numCores <- parallel::detectCores() - 1
# backend <- start_backend(numCores, cluster_type = "psock", backend_type = "async")
# 
# 
# ## load libraries in cluster
# parabar::evaluate(backend, {
#   # !all libraries that are used in the do_sim function
# })
# 
# ## load objects in cluster
# export(backend, c("cond", "do_sim"))

#### perform simulation ####
start  <- Sys.time()
test <- lapply(1:nrow(cond), do_sim, cond = cond, outputfile = "test.csv", verbose = TRUE)
#output <- par_lapply(backend, 1:nrow(cond), do_sim, cond = cond, outputfile = outputfile, verbose = FALSE)
end <- Sys.time()

# close cluster:
# stop_backend(backend)
