#### This script defines auxiliary functions for the simulation ####

#### sim_VAR() ####
# this function generates data for a single individual according to a vector autoregressive model
sim_VAR <- function(factors, obs, phi, zeta, mu, burn_in = 0){
  # obs = number of observations
  # phi = auto-regressive effect (a matrix in case of multiple constructs)
  # zeta = innovation variance (a matrix in case of multiple constructs)
  # intercept = intercept (a vector in case of multiple constructs)
  # burn_in = length of burn in (i.e., data that are generated to remove influence of initial random draw)
  
  
  # create empty dataframe of length obs + burn_in
  data <- as.data.frame(matrix(NA, nrow = burn_in + obs, ncol = factors))
  names(data) <- paste0("eta", 1:factors)
  
  intercept <- solve(solve(diag(factors) - phi), mu)
  
  for(i in 1:nrow(data)){
    # simulate the first observation from the person's starting point (= intercept)
    if(i == 1){
      data[i,] <- mu
    }
    
    # then loop through all the rows, predict the current observation from the previous observations, then add random innovation
    if(i > 1){
      predicted <- intercept + phi %*% unlist(data[i-1,])
      data[i, ] <- MASS::mvrnorm(n = 1, mu = predicted, Sigma = zeta)
    }
  }
  
  # remove the first rows, depending on length of burn in
  if(burn_in > 0){
    data <- dplyr::slice(data, burn_in+1:n()) 
  }
  
  data$obs <- 1:nrow(data)
  
  return(data)
}

#### EStep() ####
EStep <- function(pi_ks, ngroup, nclus, loglik){
  
  max_g <-rep(0,ngroup)
  z_gks <- matrix(NA,nrow = ngroup,ncol = nclus)
  
  for(g in 1:ngroup){
    for(k in 1:nclus){
      z_gks[g,k] <- log(pi_ks[k])+loglik[g,k]
    }
    max_g[g] <- max(z_gks[g,]) # prevent arithmetic underflow 
    z_gks[g,] <- exp(z_gks[g,]-rep(max_g[g],nclus))
  }
  
  # divide by the rowwise sum of the above calculated part 
  z_gks <- diag(1/apply(z_gks,1,sum))%*%z_gks
  # z_gks <- round(z_gks, digits = 16)
  # z_gks <- diag(1/apply(z_gks,1,sum))%*%z_gks
  
  return(z_gks)
}
# taken from https://github.com/AndresFPA/mmgsem/blob/main/R/E_Step.R

#### generate starting values ####
generate_startval <- function(model){
  values <- coef(model)
  values[grep("^phi", names(values))] <- runif(4, -.3, .3)
  values[c("zeta1", "zeta2")] <- runif(2, .5, 1.5)
  values["zeta12"] <- runif(1, -.3, .3)
  model <- omxSetParameters(model,
                            labels = names(values),
                            values = values)
  return(model)
}

#### create OpenMx model ####
create_model <- function(modelname, weights, objectives, model_list){
  weighted_objectives <- paste(weights, "*", objectives, collapse = " + ")
  model <- mxModel(modelname, model_list, 
                   mxAlgebraFromString(weighted_objectives, 
                                       name = "weightedfit"), 
                   mxFitFunctionAlgebra("weightedfit"))
  return(model)
}

# extract casewise log-likelihoods from fitted OpenMx models ####
get_casewiseLL <- function(models_run, n_clusters, n){
  casewiseLL <- matrix(NA,
                       nrow = n, 
                       ncol = n_clusters)
  for(i in 1:n_clusters){
    casewiseLL[, i] <- purrr::map_dbl(models_run[[i]]$submodels, ~ .x$fitfunction$result)
  }
  casewiseLL <- casewiseLL/(-2)                                                 #OpenMx gives -2 log Likelihood, so I need to divide by -2
  return(casewiseLL)
}

#### compute observed data log-likelihood
compute_observed_data_LL <- function(casewiseLL, class_proportions){
  observed_data_LL <- sum(log(rowSums(class_proportions*exp(casewiseLL))))
  return(observed_data_LL)
}

#### safely/quietly functions ####
run_step1 <- quietly(safely(step1))
run_step2 <- quietly(safely(step2))
run_step3 <- quietly(safely(step3))