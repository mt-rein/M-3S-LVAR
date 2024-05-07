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
source("step3_OpenMx_v1.R")
source("do_sim.R")
source("auxiliary_functions.R")

#### define condition grid
cond <- expand.grid(replication = 1:10,
                    n = c(48, 96),#n = c(48, 96),
                    obs = c(25, 50, 100),#obs = c(25, 50, 75, 100),
                    n_k = c(2, 4),
                    k_size = c("balanced", "unbalanced"),
                    rho_gen = c("medium", "large"),#rho_gen = c("small", "medium", "large", "very large"),
                    similarity = c("dissimilar"), #c("similar", "dissimilar")
                    innovars = "equal")#c("equal", "random")

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
})

## load objects in cluster
export(backend, c("cond", "do_sim", "run_step1", "run_step2", "run_step3",
                  "sim_VAR", "step1", "step2", "step3"))

#### perform simulation ####
start  <- Sys.time()
output <- par_lapply(backend, 1:nrow(cond), do_sim, cond = cond, outputfile = "test.csv", verbose = FALSE)
end <- Sys.time()

# close cluster:
stop_backend(backend)

test <- lapply(2:nrow(cond), do_sim, cond = cond, outputfile = "test.csv", verbose = TRUE)


library(readr)
results <- read_csv("test.csv") |> 
  arrange(iteration) |> 
  mutate(k_size = factor(k_size, levels = c("balanced", "unbalanced")),
         rho_gen = factor(rho_gen, levels = c("medium", "large")),
         similarity = factor(similarity),
         innovars = factor(innovars))

sum(results$step3_error)
