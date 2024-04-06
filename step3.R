# only for testing!
# load("examples/outputstep2.RData")
# step2output <- output
# structuralmodel = NULL

step3 <- function(step2output, structuralmodel = NULL, n_clusters,
                  nstarts = 20, maxit = 100, verbose = FALSE){
  # step2output:
  #   the object that was generated using the step2() function
  # structuralmodel (optional):
  #   a string describing the measurement model using the lavaan syntax.
  #   If omitted, a VAR model with all auto- and cross-regressive parameters
  #   will be specified by default.
  
  #### 1) Preparations ####
  ## required packages:
  # library(dplyr)
  # library(tidyr)
  # library(lavaan)
  # library(RcppAlgos)
  # library(stringr)
  
  ## extract objects from step 1 output:
  data <- step2output$data
  rho <- step2output$rho
  kappa <- step2output$kappa
  fit_step1 <- step2output$fit_step1
  
  ## generate objects for later use
  factors <- lavaan::lavNames(fit_step1, "lv")                                  # names of latent factors
  factors_ind <- paste0(factors, "_ind")                                        # names of factor score variables (single indicators)
  id <- lavaan::lavInspect(fit_step1, "cluster")                                # name of the variable that served as cluster variable in step1
  counts <- data  |>  
    dplyr::group_by(!!rlang::sym(id)) |> 
    count()
  max_t <- max(counts$n)
    
  
  #### 2) required data manipulations ####
  ## rename the factor score variables in the data
  # to use them as indicators of the latent variables
  data <- data |> 
    dplyr::rename_with(~ factors_ind, all_of(factors))
  
 if(length(factors) == 1){
    data <- data |> 
      dplyr::select(all_of(c(id, factors_ind)), k_true) |> # REMOVE K-TRUE FOR PACKAGE
      dplyr::group_by(!!rlang::sym(id)) |> 
      dplyr::mutate(obs = row_number()) |> 
      dplyr::ungroup() |> 
      tidyr::pivot_wider(names_from = obs,
                         values_from = all_of(factors_ind),
                         names_prefix = paste0(factors_ind, "_"))
  }
  if(length(factors) > 1){
    data <- data |> 
      dplyr::select(all_of(c(id, factors_ind)), k_true) |>  # REMOVE K-TRUE FOR PACKAGE
      dplyr::group_by(!!rlang::sym(id)) |> 
      dplyr::mutate(obs = row_number()) |> 
      dplyr::ungroup() |> 
      tidyr::pivot_wider(names_from = obs,
                         values_from = all_of(factors_ind))
  }
  
  #### 3) write the "measurement" model (factor measured by single indicator) ####
  # (also includes the random intercept)
  measurementmodel <- NULL
  for (fac in factors){
    indicators <- names(data)[grep(fac, names(data))]
    latents <- paste0(fac, "_", 1:length(indicators))
    
    # fix loadings to rho:
    for(col in 1:max_t){
      measurementmodel <- paste0(measurementmodel, latents[col], " =~ ", rho[fac], "*", indicators[col], " \n")
    }
    # fix residual variances to kappa:
    for(col in 1:length(indicators)){
      measurementmodel <- paste0(measurementmodel, indicators[col], " ~~ ", kappa[fac], "*", indicators[col], " \n")
    }
    
    # add population mean (fixed across time):
    measurementmodel <- paste0(measurementmodel, paste0(indicators, collapse = " + "), " ~ c(", paste0(rep(paste0("grandmean_", fac), n_clusters), collapse = ", "), ") * 1", "\n")
    # add random intercept (i.e., stable deviation of person i from grand mean):
    measurementmodel <- paste0(measurementmodel, paste0("RI_", fac), " =~ ", paste0(rho[fac], " * ", indicators, collapse = " + "), "\n")
    
    measurementmodel <- paste0(measurementmodel, paste0("RI_", fac), " ~~ ", paste0("RI_", fac), "\n")
  }
  
 #### 4) write the structural model ####
  if(is.null(structuralmodel)){
    for (fac1 in factors){
      latents1 <- paste0(fac1, "_", 1:max_t)
      for (fac2 in factors){
        latents2 <- paste0(fac2, "_", 1:max_t)
        # save 
        for(col in 2:max_t){
          # add regression effects:
          structuralmodel <- paste0(structuralmodel, latents1[col], " ~ ",
                                    paste0("c(", paste0("phi_", fac1, "_", fac2, "_k", 1:n_clusters, collapse = ", "), ") * "), latents2[col-1], " \n")
        }
      }
    }
    # add variances at first timepoint:
    
    combs <- RcppAlgos::comboGeneral(factors, 2, repetition = TRUE)
    for(comb in 1:nrow(combs)){
      for(t in 1:max_t){
        if(t == 1){
          param <- "psi_"
        } else {
          param <- "zeta_"
        }
        structuralmodel <- paste(structuralmodel,
                                 paste0(paste0(combs[comb, 1], "_", t),
                                        paste0(" ~~ c(", paste0(rep(paste0(param, combs[comb, 1], "_", combs[comb, 2]), n_clusters), collapse = ", "), ") * "),
                                        paste0(combs[comb, 2], "_", t)),
                                 "\n")
      }
    }
  }
  
  # combine models
  fullmodel <- paste(measurementmodel, structuralmodel)
  
  #### 4) preparations for clustering ####
  #### data expansion ####
  # "multiply" the data (for each cluster)
  data_expanded <- data
  for(i in 2:n_clusters){
    data_expanded <- rbind(data_expanded, data)
  }
  # create cluster assignment matrix
  cluster <- array(0, dim=c(nrow(data), n_clusters))
  for(i in 1:n_clusters){
    cluster[,i] <- i
  }
  
  cluster <- c(cluster)
  # add cluster to expanded data
  data_expanded <- cbind(data_expanded, cluster)
  
  
  #### 5) mixture modeling ####
  # maximum number of iterations:
  set.seed(8389493)#!!!! work on the whole replicability thing!!!!
  # provide seeds for the multiple starts (for replicability)
  seeds <- sample(1:100000000, nstarts)
  best_fit <- NULL
  nonconvergences <- 0
  estimation_start <- Sys.time()
  # loop across the random starts
  for(i in 1:nstarts){
    # set seed:
    set.seed(seeds[i])
    # create random posterior probabilities to start the algorithm (fake E-step)
    post <- array(runif(n_clusters*nrow(data), 0, 1), dim = c(nrow(data), n_clusters))
    post <- post/rowSums(post)
    # create weights for the data frame and append them
    w <- as.vector(post)
    data_expanded_start <- cbind(data_expanded, w)
    
    # loop over iterations until convergence or max iterations are reached
    for(it in 1:maxit) {
      # M-step
      #fitting lavaan model and updating the class proportions
      if(it == 1) {
        output <- lavaan(fullmodel, data = data_expanded_start, int.ov.free = FALSE, auto.var = FALSE,
                             group = "cluster", sampling.weight = "w",
                             baseline = FALSE, se = "none", h1 = FALSE)
        fit <- output
        start <- coef(fit)
      } else {
        output <- lavaan(fullmodel, data = data_expanded_start, int.ov.free = FALSE, auto.var = FALSE,
                             group = "cluster", sampling.weight = "w",
                             baseline = FALSE, se = "none", h1 = FALSE,
                             check.post = FALSE,
                             control = list(rel.tol = 1e-06),
                             start = start)
        fit <- output
        start <- coef(fit)
      }
      # compute class proportions:
      pi = colMeans(post)
      
      #E-step
      #computing new posteriors and replacing these in the expanded data file
      loglik.case <- unlist(lavInspect(fit, what = "loglik.casewise"), use.names = FALSE) / (w * n_clusters)
      # note: the casewise logliks are scaled by number of groups (i.e., clusters) and weight, hence they are "unweighed" here
      loglik.case <- matrix(loglik.case, ncol = n_clusters)
      post <- t(t(exp(loglik.case))*pi)
      lik <- rowSums(post)
      post <- post/lik
      
      
      # from posteriors, get new weights:
      w <- as.vector(post)
      data_expanded_start$w <- w
      
      #calculating log-likelihood value and checking convergence
      loglik=sum(log(lik))
      if(verbose){
        print(c(it,loglik))
        }
      if(it>1&&(loglik-loglik0)<1.0e-6){
        break
      }
      loglik0=loglik
      if(it == maxit){
        if(verbose){
          print(paste("Max iterations reached without convergence. Start:", i))
        }
        nonconvergences <- nonconvergences + 1
      }
    }
    
    # check if the new fit is better than the best one from previous random starts
    if(i == 1){
      # if it's the first random start, save results
      best_loglik <- loglik
      best_w <- w
    } else {
      # otherwise, compare fits
      if(loglik > best_loglik){
        best_loglik <- loglik
        best_w <- w
      }
    }
    
    if(verbose){
      print(paste("Start", i, "out of", nstarts, "completed."))
      }
  }
  duration <- difftime(Sys.time(), estimation_start, unit = "s")
  
  # fit one more time to get standard errors
  w <- best_w
  data_expanded <- cbind(data_expanded, w)
  fit <- lavaan(fullmodel, data = data_expanded, int.ov.free = FALSE, auto.var = FALSE,
                group = "cluster", sampling.weight = "w", ridge = TRUE)
  
  #### find proxy maximum ###
  # w_true <- ifelse(data_expanded$cluster == data_expanded$k_true, .9999999, 1-.9999999) # this doesn't work! Need less extreme values I think, then feed to algorithm?
  # data_expanded$w <- w_true
  # post_true <- matrix(w_true, ncol = n_clusters)
  # fit_true <- lavaan(fullmodel, data = data_expanded, int.ov.free = FALSE, auto.var = FALSE,
  #                    group = "cluster", sampling.weight = "w", ridge = TRUE)

  
  #### 5) extract estimates ####
  ## TO DO: make this more flexible to accommodate for custom SM (provided by user)\
  params <- coef(fit)
  phi_names <- unique(names(params))[grep("^phi", unique(names(params)))]
  phis <- purrr::map_dbl(phi_names,\(x) params[names(params) == x][1])
  names(phis) <- phi_names
  
  psi_names <- unique(names(params))[grep("^psi", unique(names(params)))]
  psis <- purrr::map_dbl(psi_names,\(x) params[names(params) == x][1])
  names(psis) <- psi_names
  
  zeta_names <- unique(names(params))[grep("^zeta", unique(names(params)))]
  zetas <- purrr::map_dbl(zeta_names,\(x) params[names(params) == x][1])
  names(zetas) <- zeta_names
  
  estimates <- list("phi" = phis,
                    "psi" = psis,
                    "zeta" = zetas)
  
  loglik <- best_loglik
  post <- matrix(best_w, ncol = 2)
  colnames(post) <- 1:n_clusters
  assignment <- apply(post, MARGIN = 1, function(x) ifelse(x == max(x), 1, 0))
  assignment <- t(assignment)
  pi <- colMeans(post)
  
  clustering <- list("pi" = pi,
                     "posterior_prob" = post,
                     "assignment" = assignment)
  
  other <- list("loglik" = loglik,
                "nonconvergences" = nonconvergences,
                #"proxy_maximum" = proxy_maximum, #remove for package
                "duration" = duration)
  
  #### 5) build the output ####
  output <- list("fit_step3" = fit,
                 "data" = data,
                 "estimates" = estimates,
                 "clustering" = clustering,
                 "other" = other)
  
  return(output)
}