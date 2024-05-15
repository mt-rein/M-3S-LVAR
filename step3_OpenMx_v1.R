# only for testing!
# load("examples/outputstep2.RData")
# step2output <- output
# structuralmodel = NULL

step3 <- function(step2output, n_clusters, n_starts = 25, n_best_starts = 5,
                  maxit = 100,
                  true_clusters, verbose = FALSE){
  # step2output:
  #   the object that was generated using the step2() function
  # ...
  
  
  #### FOR TESTING ####
  # step2output = output_step2$result$result
  # n_clusters = n_k
  # n_starts = 25
  # n_best_starts = 5
  # maxit = 100
  # true_clusters = clusterassignment_true
  # verbose = TRUE
  
  #### 1) Preparations ####
  ## extract objects from step 1 output:
  data <- step2output$data
  rho <- step2output$rho
  kappa <- step2output$kappa
  fit_step1 <- step2output$fit_step1
  
  ## generate objects for later use
  factors <- lavaan::lavNames(fit_step1, "lv")                                  # names of latent factors
  factors_ind <- paste0(factors, "_ind")                                        # names of factor score variables (single indicators)
  id <- lavaan::lavInspect(fit_step1, "cluster")                                # name of the variable that served as cluster variable in step1
  unique_ids <- unique(data[, id])                                              # vector of unique ids
  n <- length(unique_ids)                                                       # number of individuals
  
  # save previous options of OpenMx, then set new ones:
  mx_Hessian <- mxOption(key = "Calculate Hessian")
  mx_SE <- mxOption(key = "Standard Errors")
  mx_maxiterations <- mxOption(key = "Major Iterations")
  mxOption(key = "Calculate Hessian", value = "No")
  mxOption(key = "Standard Errors", value = "No")
  
  
  #### 2) data manipulation ####
  ## rename the factor score variables in the data
  # to use them as indicators of the latent variables
  data <- data |> 
    dplyr::rename_with(~ factors_ind, all_of(factors))
  
  #### 3) create OpenMx matrices ####
  xdim <- length(factors)*2 # number of latent constructs in model is number of factors times 2 due to the random intercepts
  udim <- 1 # exogenous covariates (ignored so far)
  ydim <- length(factors) # number of indicators (i.e., factor score variables)
  
  ## A matrices (= dynamics)
  amat <- mxMatrix(type = "Full", nrow = xdim, ncol = xdim,
                   free = c(TRUE, TRUE, FALSE, FALSE,
                            TRUE, TRUE, FALSE, FALSE,
                            FALSE, FALSE, FALSE, FALSE,
                            FALSE, FALSE, FALSE, FALSE),
                   values = c(.1, .1, 1, 0,
                              .1, .1, 0, 1,
                              0, 0, 1, 0,
                              0, 0, 0, 1),
                   name = "A",
                   labels = c("phi11", "phi12", NA, NA,
                              "phi21", "phi22", NA, NA,
                              NA, NA, NA, NA,
                              NA, NA, NA, NA),
                   lbound= c(-.9, -.9, NA, NA,
                             -.9, -.9, NA, NA,
                             NA, NA, NA, NA,
                             NA, NA, NA, NA),
                   ubound= c(.9, .9, NA, NA,
                             .9, .9, NA, NA,
                             NA, NA, NA, NA,
                             NA, NA, NA, NA),
                   byrow = TRUE)
  
  # B matrix (= exogenous covariates on latent constructs)
  bmat <- mxMatrix('Zero', nrow = xdim, ncol = udim,
                   name='B')
  
  # C matrix (= factor loadings, here fixed to rho)
  cmat <- mxMatrix('Full', nrow = ydim, ncol = xdim,
                   free = c(FALSE),
                   values = c(rho[1], 0, 0, 0,
                              0, rho[2], 0, 0),
                   labels = NA,
                   byrow = TRUE,
                   name='C',
                   dimnames = list(factors_ind, c(paste0(factors),
                                                  paste0("intercept_", factors))
                                   )
                   )
  
  # D matrix (= exogenous covariates on observed variables)
  dmat <- mxMatrix('Zero', nrow = ydim, ncol = udim,
                   name='D')
  
  # Q matrix (= innovation (co)variances)
  qmat <- mxMatrix("Full", xdim, xdim,
                   free = c(TRUE, TRUE, FALSE, FALSE,
                            TRUE, TRUE, FALSE, FALSE,
                            FALSE, FALSE, FALSE, FALSE,
                            FALSE, FALSE, FALSE, FALSE),
                   values = c(1, .3, 0, 0,
                              .3, 1, 0, 0,
                              0, 0,  0, 0,
                              0, 0,  0, 0),
                   name = "Q",
                   labels = c("zeta1", "zeta12", NA, NA,
                              "zeta12", "zeta2", NA, NA,
                              NA, NA, NA, NA,
                              NA, NA, NA, NA),
                   lbound= c(1e-6, NA, NA, NA,
                             NA, 1e-6, NA, NA,
                             NA, NA, NA, NA,
                             NA, NA, NA, NA),
                   byrow = TRUE)
  
  
  # R matrix (= measurement noise, here fixed to kappa)
  rmat <- mxMatrix('Diag', nrow = ydim, ncol = ydim,
                   free = FALSE,
                   values = kappa,
                   name = 'R')
  
  # x0 and P0 (= initial values and (co)variances of the latent constructs)
  xmat <- mxMatrix('Full', nrow = xdim, ncol = 1,
                   free = TRUE,
                   values = c(0, 0, 1, 1),
                   name='x0',
                   labels = c(paste0("ini_", factors),
                              paste0("m_intercept_", factors))
                   )
  
  pmat <- mxMatrix('Diag', nrow = xdim, ncol = xdim,
                   free = TRUE,
                   values = c(1, 1, 1, 1),
                   name='P0',
                   labels = c(paste0("var_", factors),
                              paste0("var_intercept_", factors))
                   )
  
  # u (= covariates)
  umat <- mxMatrix('Zero', nrow = udim, ncol = 1, name='u')
  
  #### 4) create OpenMx models ####
  # create a list of models (one for each individual) for each latent class:
  modelnames <- paste0("id_", unique_ids)
  objectives <- paste0(modelnames, ".objective")                                # vector of objective strings, needed later to create custom fit functions
  fitfunctions <- paste0(modelnames, ".fitfunction")                            # vector of fit function strings, needed later to extract case wise log-likelihoods
  model_list_k1 <- vector(mode = "list", length = n)
  for(i in unique_ids){
    model_list_k1[[i]] <- mxModel(name = modelnames[i],
                                  amat, bmat, cmat, dmat,
                                  qmat, rmat, xmat, pmat, 
                                  umat,
                                  mxExpectationStateSpace('A', 'B', 'C', 'D', 
                                                          'Q', 'R', 'x0', 'P0', 
                                                          'u'),
                                  mxFitFunctionML(),
                                  mxData(data[data[, id] == i, factors_ind], 
                                         'raw'))
  }
  
  model_list_k2 <- vector(mode = "list", length = n)
  for(i in unique_ids){
    model_list_k2[[i]] <- mxModel(name = modelnames[i],
                                  amat, bmat, cmat, dmat, 
                                  qmat, rmat, xmat, pmat, 
                                  umat,
                                  mxExpectationStateSpace('A', 'B', 'C', 'D', 
                                                          'Q', 'R', 'x0', 'P0', 
                                                          'u'),
                                  mxFitFunctionML(),
                                  mxData(data[data[, id] == i, factors_ind], 
                                         'raw'))
    }
  
  if(n_clusters == 4){
    model_list_k3 <- vector(mode = "list", length = n)
    for(i in unique_ids){
      model_list_k3[[i]] <- mxModel(name = modelnames[i],
                                    amat, bmat, cmat, dmat, 
                                    qmat, rmat, xmat, pmat, 
                                    umat,
                                    mxExpectationStateSpace('A', 'B', 'C', 'D', 
                                                            'Q', 'R', 'x0', 'P0', 
                                                            'u'),
                                    mxFitFunctionML(),
                                    mxData(data[data[, id] == i, factors_ind], 
                                           'raw'))
    }
    
    model_list_k4 <- vector(mode = "list", length = n)
    for(i in unique_ids){
      model_list_k4[[i]] <- mxModel(name = modelnames[i],
                                    amat, bmat, cmat, dmat, qmat, rmat, xmat, pmat, umat,
                                    mxExpectationStateSpace('A', 'B', 'C', 'D', 
                                                            'Q', 'R', 'x0', 'P0', 
                                                            'u'),
                                    mxFitFunctionML(),
                                    mxData(data[data[, id] == i, factors_ind], 
                                           'raw'))
      }
  }
  
  #### 5) mixture modeling ####
  #set.seed(8389493)#!!!! work on the whole replicability thing!!!!
  # provide seeds for the multiple starts (for replicability)
  seeds <- sample(1:100000000, n_starts)
  estimation_start <- Sys.time()
  # The estimation happens in two parts: First, a single iteration of the EM algorithm is run for each random start
  # the n_best_starts best starts are then completed
  
  ## Part 1: all random starts
  all_starts <- vector(mode = "list", length = n_starts)                        # create empty vector for the output
  mxOption(key = "Major iterations", value = 3)                                 # reduce the number of iterations when estimating the parameters in OpenMx (to speed it up)
  for(i in 1:n_starts){
    # set seed:
    set.seed(seeds[i])
    # create random cluster assignment:
    post <- matrix(0, nrow = n, ncol = n_clusters)
    for(person in 1:n){
      if(person <= n_clusters){
        # assign the first persons to certain clusters to avoid empty clusters
        # (i.e., first person in first cluster, second person in second cluster, etc)
        # until every cluster has one individual
        post[person, person] <- 1
      } else {
        post[person, sample(1:n_clusters, 1)] <- 1
      }
    }
    
    # compute class proportions:
    class_proportions <- colMeans(post)
    
    # fitting SSM model and update parameter estimates
    ## cluster 1:
    # create custom fit function that corresponds to complete data log-likelihood:
    weights = post[, 1]
    weighted_objectives <- paste(weights, "*", objectives, collapse = " + ")
    # create model from person-models:
    model_k1 <- mxModel('model_k1', model_list_k1, 
                        mxAlgebraFromString(weighted_objectives, 
                                            name = "weightedfit"), 
                        mxFitFunctionAlgebra("weightedfit"))
    # generate random starting values: 
    model_k1 <- generate_startval(model_k1)
      
    # run the model:
    model_k1r <- mxRun(model_k1, silent = !verbose, suppressWarnings = TRUE)
    
    ## cluster 2:
    weights = post[, 2]
    weighted_objectives <- paste(weights, "*", objectives, collapse = " + ")
    model_k2 <- mxModel('model_k2', model_list_k2, 
                        mxAlgebraFromString(weighted_objectives, 
                                            name = "weightedfit"), 
                        mxFitFunctionAlgebra("weightedfit"))
    model_k2 <- generate_startval(model_k2)
    model_k2r <- mxRun(model_k2, silent = !verbose, suppressWarnings = TRUE)
    
    # create models for clusters 3 and 4 (if applicable):
    if(n_clusters == 4){
      ## cluster 3:
      weights = post[, 3]
      weighted_objectives <- paste(weights, "*", objectives, collapse = " + ")
      model_k3 <- mxModel('model_k3', model_list_k3, 
                          mxAlgebraFromString(weighted_objectives, 
                                              name = "weightedfit"), 
                          mxFitFunctionAlgebra("weightedfit"))
      model_k3 <- generate_startval(model_k3)
      model_k3r <- mxRun(model_k3, silent = !verbose, suppressWarnings = TRUE)
      
      ## cluster 4:
      weights = post[, 4]
      weighted_objectives <- paste(weights, "*", objectives, collapse = " + ")
      model_k4 <- mxModel('model_k4', model_list_k4, 
                          mxAlgebraFromString(weighted_objectives, 
                                              name = "weightedfit"), 
                          mxFitFunctionAlgebra("weightedfit"))
      model_k4 <- generate_startval(model_k4)
      model_k4r <- mxRun(model_k4, silent = !verbose, suppressWarnings = TRUE)
    }
    
    ## compute LL:
    m2_loglik_k1 <- purrr::map_dbl(fitfunctions, mxEvalByName, model=model_k1r)
    m2_loglik_k2 <- purrr::map_dbl(fitfunctions, mxEvalByName, model=model_k2r)
    if(n_clusters == 4){
      m2_loglik_k3 <- purrr::map_dbl(fitfunctions, mxEvalByName, model=model_k3r)
      m2_loglik_k4 <- purrr::map_dbl(fitfunctions, mxEvalByName, model=model_k4r)
    }
    # combine into a matrix with n rows and k columns:
    if(n_clusters == 2){
      casewise.m2_loglik <- cbind(m2_loglik_k1, m2_loglik_k2) 
    }
    if(n_clusters == 4){
      casewise.m2_loglik <- cbind(m2_loglik_k1, m2_loglik_k2, m2_loglik_k3, m2_loglik_k4) 
    }
    # transform on log likelihood scale:
    casewise.loglik <- casewise.m2_loglik/(-2)
    
    # compute complete-data log likelihood:
    # complete_data_loglik <-  sum(post*log(class_proportions)-post*casewise.loglik)
    
    # compute observed-data log likelihood:
    observed_data_loglik <- sum(log(rowSums(class_proportions*exp(casewise.loglik))))
    
    # save observed_data_loglik, post, and parameter estimates to use as starting values in Part 2:
    if(n_clusters == 2){
      all_starts[[i]] <- list("loglik" = observed_data_loglik,
                              "post" = post,
                              "coefs" = list("cluster1" = coef(model_k1),
                                             "cluster2" = coef(model_k2)
                                             )
                              )
    }
    if(n_clusters == 4){
      all_starts[[i]] <- list("loglik" = observed_data_loglik,
                              "post" = post,
                              "coefs" = list("cluster1" = coef(model_k1),
                                             "cluster2" = coef(model_k2),
                                             "cluster3" = coef(model_k3),
                                             "cluster4" = coef(model_k4)
                                             )
                              )
    }
  }
  
  if(verbose){
    print("Setup finished. Finishing 5 best starts now.")
  }
  
  ## Part 2: complete the best starts
  # extract the best starts:
  best_starts <- purrr::map_dbl(all_starts, ~ .x$loglik) |> 
    order(decreasing = TRUE) |> 
    head(n_best_starts)
  
  best_loglik <- -Inf
  most_iterations <- 0
  nonconvergences <- 0
  for(i in best_starts){
    mxOption(key = "Major iterations", value = 3)                               # reset the maximum number of iterations in OpenMx to 3
    close_convergence <- FALSE                                                  # reset the close_convergence switch
    
    # extract cluster assignment from corresponding start in Part 1:
    post <- all_starts[[i]]$post
    
    # loop over iterations until convergence or max iterations are reached
    for(it in 1:maxit){
      ## E-step: update class membership (skip in first iteration) and class proportions
      if(it >1){
        # update posteriors:
        post <- EStep(pi = class_proportions, ngroup = n, 
                      nclus = n_clusters, loglik = casewise.loglik)
      }
      
      # compute class proportions:
      class_proportions <- colMeans(post)
      
      ## M-step: fitting SSM model and update parameter estimates
      # same procedure as above for each cluster
      # cluster 1:
      weights = post[, 1]
      weighted_objectives <- paste(weights, "*", objectives, collapse = " + ")
      model_k1 <- mxModel('model_k1', model_list_k1, 
                          mxAlgebraFromString(weighted_objectives, 
                                              name = "weightedfit"), 
                          mxFitFunctionAlgebra("weightedfit")
                          )
      if(it == 1){
        # in first iteration, get starting values from corresponding start in Part 1:
        model_k1 <- omxSetParameters(model_k1, 
                                     labels = names(all_starts[[i]]$coefs$cluster1), 
                                     values = all_starts[[i]]$coefs$cluster1)
      }
      if(it > 1){
        # in following iterations, extract starting values from estimates of previous iteration:
        model_k1 <- omxSetParameters(model_k1, 
                                     labels = names(coef(model_k1r)), 
                                     values = coef(model_k1r))
      }
      model_k1r <- mxRun(model_k1, silent = !verbose, suppressWarnings = TRUE)
      
      # cluster 2:
      weights = post[, 2]
      weighted_objectives <- paste(weights, "*", objectives, collapse = " + ")
      model_k2 <- mxModel('model_k2', model_list_k2, 
                          mxAlgebraFromString(weighted_objectives, 
                                              name = "weightedfit"), 
                          mxFitFunctionAlgebra("weightedfit")
                          )
      if(it == 1){
        model_k2 <- omxSetParameters(model_k2, 
                                     labels = names(all_starts[[i]]$coefs$cluster2), 
                                     values = all_starts[[i]]$coefs$cluster2)
      }
      if(it > 1){
        model_k2 <- omxSetParameters(model_k2, 
                                     labels = names(coef(model_k2r)), 
                                     values = coef(model_k2r))
      }
      model_k2r <- mxRun(model_k2, silent = !verbose, suppressWarnings = TRUE)
      
      # create models for clusters 3 and 4 (if applicable):
      if(n_clusters == 4){
        # cluster 3:
        weights = post[, 3]
        weighted_objectives <- paste(weights, "*", objectives, collapse = " + ")
        model_k3 <- mxModel('model_k3', model_list_k3, 
                            mxAlgebraFromString(weighted_objectives, 
                                                name = "weightedfit"), 
                            mxFitFunctionAlgebra("weightedfit"))
        if(it == 1){
          model_k3 <- omxSetParameters(model_k3, 
                                       labels = names(all_starts[[i]]$coefs$cluster3), 
                                       values = all_starts[[i]]$coefs$cluster3)
        }
        if(it > 1){
          model_k3 <- omxSetParameters(model_k3, 
                                       labels = names(coef(model_k3r)), 
                                       values = coef(model_k3r))
        }
        model_k3r <- mxRun(model_k3, silent = !verbose, suppressWarnings = TRUE)
        
        # cluster 4:
        weights = post[, 4]
        weighted_objectives <- paste(weights, "*", objectives, collapse = " + ")
        model_k4 <- mxModel('model_k4', model_list_k4, 
                            mxAlgebraFromString(weighted_objectives, 
                                                name = "weightedfit"), 
                            mxFitFunctionAlgebra("weightedfit"))
        if(it == 1){
          model_k4 <- omxSetParameters(model_k4, 
                                       labels = names(all_starts[[i]]$coefs$cluster4), 
                                       values = all_starts[[i]]$coefs$cluster4)
        }
        if(it > 1){
          model_k4 <- omxSetParameters(model_k4, 
                                       labels = names(coef(model_k4r)), 
                                       values = coef(model_k4r))
        }
        model_k4r <- mxRun(model_k4, silent = !verbose, suppressWarnings = TRUE)
      }
      
      ## compute new LL:
      fitfunctions <- paste0(modelnames, ".fitfunction")
      m2_loglik_k1 <- purrr::map_dbl(fitfunctions, mxEvalByName, model=model_k1r)
      m2_loglik_k2 <- purrr::map_dbl(fitfunctions, mxEvalByName, model=model_k2r)
      if(n_clusters == 4){
        m2_loglik_k3 <- purrr::map_dbl(fitfunctions, mxEvalByName, model=model_k3r)
        m2_loglik_k4 <- purrr::map_dbl(fitfunctions, mxEvalByName, model=model_k4r)
      }
      # combine into a matrix with n rows and k columns
      if(n_clusters == 2){
        casewise.m2_loglik <- cbind(m2_loglik_k1, m2_loglik_k2) 
      }
      if(n_clusters == 4){
        casewise.m2_loglik <- cbind(m2_loglik_k1, m2_loglik_k2, m2_loglik_k3, m2_loglik_k4) 
      }
      # transform on log likelihood scale:
      casewise.loglik <- casewise.m2_loglik/(-2)
      
      
      # compute complete-data log likelihood
      # complete_data_loglik <-  sum(post*log(class_proportions)-post*casewise.loglik)
      
      # compute observed-data log likelihood
      observed_data_loglik <- sum(log(rowSums(class_proportions*exp(casewise.loglik))))
      
      if(verbose){
        if(it != 1){
          print(paste0("Iteration: ", it, ". Log Likelihood: ", round(observed_data_loglik, 4), ". Change: ", round(observed_data_loglik - observed_data_loglik0, 6), "."))
        } else {
          print(paste0("Iteration: ", it, ". Log Likelihood: ", round(observed_data_loglik, 4), "."))
        }
        
      }
      
      # check convergence and break loop if applicable:
      if(it > 2 && abs(observed_data_loglik - observed_data_loglik0) < 1.0e-6){
        if(verbose){
          print("Convergence achieved.")
        }
        # save models to extract estimates:
        if(n_clusters == 2){
          models <- list(model_k1r, model_k2r)
        }
        if(n_clusters == 4){
          models <- list(model_k1r, model_k2r, model_k3r, model_k4r)
        }
        
        break
      }
      
      # check close convergence to increase number of iterations in OpenMx:
      if(it > 5 && (observed_data_loglik - observed_data_loglik0) < 5 && !close_convergence){
        mxOption(key = "Major iterations", value = mx_maxiterations)
        close_convergence <- TRUE
        print("Close to convergence.")
      }
      
      observed_data_loglik0 = observed_data_loglik
      
      # check if maximum umber of iterations has been reached:
      if(it == maxit){
        if(verbose){
          print(paste("Max iterations reached without convergence. Start:", i))
        }
        nonconvergences <- nonconvergences + 1
      }
    }
    
    # check if the new fit is better than the best one from previous starts
    if(observed_data_loglik > best_loglik){
      best_loglik <- observed_data_loglik
      best_post <- post
      best_models <- models
    }
    if(it > most_iterations){
      most_iterations <- it
    }
    
    if(verbose){
      print(paste("Start", i, "out of", n_starts, "completed."))
    }
  }
  duration <- difftime(Sys.time(), estimation_start, unit = "s")
  
  #### find proxy maximum ###
  # create matrix with true cluster assignments
  post <- matrix(0, nrow = n, ncol = n_clusters)
  for (i in 1:n_clusters) {
    post[true_clusters == i, i] <- 1
  }
  mxOption(key = "Major iterations", value = mx_maxiterations)
  
  for(it in 1:maxit) {
    ## E-step: update class membership (skip in first iteration) and class proportions
    if(it >1){
      # update posteriors:
      post <- EStep(pi = class_proportions, ngroup = n, 
                    nclus = n_clusters, loglik = casewise.loglik)
    }
    
    # compute class proportions:
    class_proportions <- colMeans(post)
    
    ## M-step: fitting SSM model and update parameter estimates
    # cluster 1:
    weights = post[, 1]
    weighted_objectives <- paste(weights, "*", objectives, collapse = " + ")
    model_k1 <- mxModel('model_k1', model_list_k1, 
                        mxAlgebraFromString(weighted_objectives, 
                                            name = "weightedfit"), 
                        mxFitFunctionAlgebra("weightedfit"))
    if(it == 1){
      # in first iteration, get starting values from best model 
      model_k1 <- omxSetParameters(model_k1, 
                                   labels = names(coef(model_k1)), 
                                   values = coef(best_models[[1]]))
    }
    if(it > 1){
      # in following iterations, get starting values from previous iteration:
      model_k1 <- omxSetParameters(model_k1, 
                                   labels = names(coef(model_k1)), 
                                   values = coef(model_k1r))
    }
    model_k1r <- mxRun(model_k1, silent = !verbose, suppressWarnings = TRUE)
    
    # cluster 2:
    weights = post[, 2]
    weighted_objectives <- paste(weights, "*", objectives, collapse = " + ")
    model_k2 <- mxModel('model_k2', model_list_k2, 
                        mxAlgebraFromString(weighted_objectives, 
                                            name = "weightedfit"), 
                        mxFitFunctionAlgebra("weightedfit"))
    if(it == 1){
      model_k2 <- omxSetParameters(model_k2, 
                                   labels = names(coef(model_k2)), 
                                   values = coef(best_models[[2]]))
    }
    if(it > 1){
      model_k2 <- omxSetParameters(model_k2, 
                                   labels = names(coef(model_k2)), 
                                   values = coef(model_k2r))
    }
    model_k2r <- mxRun(model_k2, silent = !verbose, suppressWarnings = TRUE)
    
    # create models for clusters 3 and 4 (if applicable):
    if(n_clusters == 4){
      # cluster 3:
      weights = post[, 3]
      weighted_objectives <- paste(weights, "*", objectives, collapse = " + ")
      model_k3 <- mxModel('model_k3', model_list_k3, 
                          mxAlgebraFromString(weighted_objectives, 
                                              name = "weightedfit"), 
                          mxFitFunctionAlgebra("weightedfit"))
      if(it == 1){
        model_k3 <- omxSetParameters(model_k3, 
                                     labels = names(coef(model_k3)), 
                                     values = coef(best_models[[3]]))
      }
      if(it > 1){
        model_k3 <- omxSetParameters(model_k3, labels = names(coef(model_k3)), values = coef(model_k3r))
      }
      model_k3r <- mxRun(model_k3, silent = !verbose, suppressWarnings = TRUE)
      
      # cluster 4:
      weights = post[, 4]
      weighted_objectives <- paste(weights, "*", objectives, collapse = " + ")
      model_k4 <- mxModel('model_k4', model_list_k4, 
                          mxAlgebraFromString(weighted_objectives, 
                                              name = "weightedfit"), 
                          mxFitFunctionAlgebra("weightedfit"))
      if(it == 1){
        model_k4 <- omxSetParameters(model_k4, 
                                     labels = names(coef(model_k4)), 
                                     values = coef(best_models[[4]]))
      }
      if(it > 1){
        model_k4 <- omxSetParameters(model_k4, 
                                     labels = names(coef(model_k4)), 
                                     values = coef(model_k4r))
      }
      model_k4r <- mxRun(model_k4, silent = !verbose, suppressWarnings = TRUE)
    }
    
    ## update LL:
    m2_loglik_k1 <- purrr::map_dbl(fitfunctions, mxEvalByName, model=model_k1r)
    m2_loglik_k2 <- purrr::map_dbl(fitfunctions, mxEvalByName, model=model_k2r)
    if(n_clusters == 4){
      m2_loglik_k3 <- purrr::map_dbl(fitfunctions, mxEvalByName, model=model_k3r)
      m2_loglik_k4 <- purrr::map_dbl(fitfunctions, mxEvalByName, model=model_k4r)
    }
    # combine them into a matrix with n rows and k columns
    if(n_clusters == 2){
      casewise.m2_loglik <- cbind(m2_loglik_k1, m2_loglik_k2) 
    }
    if(n_clusters == 4){
      casewise.m2_loglik <- cbind(m2_loglik_k1, m2_loglik_k2, m2_loglik_k3, m2_loglik_k4) 
    }
    # transform on log likelihood scale:
    casewise.loglik <- casewise.m2_loglik/(-2)
    
    
    # compute complete-data log likelihood
    # complete_data_loglik <-  sum(post*log(class_proportions)-post*casewise.loglik)
    
    # compute observed-data log likelihood
    observed_data_loglik <- sum(log(rowSums(class_proportions*exp(casewise.loglik))))
    
    if(verbose){
      if(it != 1){
        print(paste0("Iteration: ", it, ". Log Likelihood: ", round(observed_data_loglik, 4), ". Change: ", round(observed_data_loglik - observed_data_loglik0, 6), "."))
      } else {
        print(paste0("Iteration: ", it, ". Log Likelihood: ", round(observed_data_loglik, 4), "."))
      }
      
    }
    
    
    # check convergence and break loop if applicable:
    if(it > 1 && (observed_data_loglik - observed_data_loglik0) < 1.0e-6){
      if(verbose){
        print("Convergence achieved.")
      }
      if(n_clusters == 2){
        models <- list(model_k1r, model_k2r)
      }
      if(n_clusters == 4){
        models <- list(model_k1r, model_k2r, model_k3r, model_k4r)
      }
      
      break
    }
    
    observed_data_loglik0 = observed_data_loglik
  }
  proxy_maximum <- observed_data_loglik
  
  
  #### 5) extract estimates ####
  estimates <- lapply(best_models, coef)
  loglik <- best_loglik
  post <- best_post
  names(estimates) <- colnames(post) <- paste0("cluster", 1:n_clusters)
  assignment <- round(post)
  class_proportions <- colMeans(post)
  
  
  
  clustering <- list("class_proportions" = class_proportions,
                     "posterior_prob" = post,
                     "assignment" = assignment)
  
  other <- list("loglik" = loglik,
                "nonconvergences" = nonconvergences,
                "proxy_maximum" = proxy_maximum,
                "duration" = duration,
                "most_iterations" = most_iterations)
  
  #### 5) build the output ####
  output <- list("data" = data,
                 "estimates" = estimates,
                 "clustering" = clustering,
                 "other" = other)
  
  ## reset the changed OpenMx Options:
  mxOption(key = "Calculate Hessian", value = mx_Hessian) 
  mxOption(key = "Standard Errors", value = mx_SE)
  mxOption(key = "Major Iterations", value = mx_maxiterations)
  
  return(output)
}