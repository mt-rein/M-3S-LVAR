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
  
  
  #### FOR TESTING ####
  # step2output = output_step2$result$result
  # n_clusters = n_k
  # nstarts = 1
  # maxit = 100
  # structuralmodel = NULL
  # verbose = TRUE
  # # 
  
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
  
  
  #### 2) data manipulation ####
  ## rename the factor score variables in the data
  # to use them as indicators of the latent variables
  data <- data |> 
    dplyr::rename_with(~ factors_ind, all_of(factors))
  
  #### 3) create OpenMx matrices ####
  xdim <- length(factors)*2 # number of latent constructs in model is number of factors times 2 due to the random intercepts
  udim <- 1 # exogenous covariates (ignored so far)
  ydim <- length(factors) # number of indicators (i.e., factor score variables)
  
  sv_k1 <- runif(4, -.3, .3)
  sv_k2 <- runif(4, -.3, .3)
  if(n_clusters == 4){
    sv_k3 <- runif(4, -.3, .3)
    sv_k4 <- runif(4, -.3, .3)
  }
  
  ## A matrices (= dynamics)
  amat_k1 <- mxMatrix(type = "Full", nrow = xdim, ncol = xdim,
                      free = c(TRUE, TRUE, FALSE, FALSE,
                               TRUE, TRUE, FALSE, FALSE,
                               FALSE, FALSE, FALSE, FALSE,
                               FALSE, FALSE, FALSE, FALSE),
                      values = c(sv_k1[1], sv_k1[2], 1, 0,
                                 sv_k1[3], sv_k1[4], 0, 1,
                                 0, 0, 1, 0,
                                 0, 0, 0, 1),
                      name = "A1",
                      labels = c("phi11_k1", "phi12_k1", NA, NA,
                                 "phi21_k1", "phi22_k1", NA, NA,
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
  
  amat_k2 <- mxMatrix(type = "Full", nrow = xdim, ncol = xdim,
                      free = c(TRUE, TRUE, FALSE, FALSE,
                               TRUE, TRUE, FALSE, FALSE,
                               FALSE, FALSE, FALSE, FALSE,
                               FALSE, FALSE, FALSE, FALSE),
                      values = c(sv_k2[1], sv_k2[2], 1, 0,
                                 sv_k2[3], sv_k2[4], 0, 1,
                                 0, 0, 1, 0,
                                 0, 0, 0, 1),
                      name = "A2",
                      labels = c("phi11_k2", "phi12_k2", NA, NA,
                                 "phi21_k2", "phi22_k2", NA, NA,
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
  
  if(n_clusters == 4){
    amat_k3 <- mxMatrix(type = "Full", nrow = xdim, ncol = xdim,
                        free = c(TRUE, TRUE, FALSE, FALSE,
                                 TRUE, TRUE, FALSE, FALSE,
                                 FALSE, FALSE, FALSE, FALSE,
                                 FALSE, FALSE, FALSE, FALSE),
                        values = c(sv_k3[1], sv_k3[2], 1, 0,
                                   sv_k3[3], sv_k3[4], 0, 1,
                                   0, 0, 1, 0,
                                   0, 0, 0, 1),
                        name = "A3",
                        labels = c("phi11_k3", "phi12_k3", NA, NA,
                                   "phi21_k3", "phi22_k3", NA, NA,
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
    
    amat_k4 <- mxMatrix(type = "Full", nrow = xdim, ncol = xdim,
                        free = c(TRUE, TRUE, FALSE, FALSE,
                                 TRUE, TRUE, FALSE, FALSE,
                                 FALSE, FALSE, FALSE, FALSE,
                                 FALSE, FALSE, FALSE, FALSE),
                        values = c(sv_k4[1], sv_k4[2], 1, 0,
                                   sv_k4[3], sv_k4[4], 0, 1,
                                   0, 0, 1, 0,
                                   0, 0, 0, 1),
                        name = "A4",
                        labels = c("phi11_k4", "phi12_k4", NA, NA,
                                   "phi21_k4", "phi22_k4", NA, NA,
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
  }
  
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
                   dimnames = list(factors_ind, c(paste0(factors), paste0("intercept_", factors))))
  
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
                   labels = c(paste0("ini_", factors), paste0("m_intercept_", factors)))
  
  pmat <- mxMatrix('Diag', nrow = xdim, ncol = xdim,
                   free = TRUE,
                   values = c(1, 1, 1, 1),
                   name='P0',
                   labels = c(paste0("var_", factors), paste0("var_intercept_", factors)))
  
  # u (= covariates)
  umat <- mxMatrix('Zero', nrow = udim, ncol = 1, name='u')
  
  # collect all invariant matrices (i.e., all except the A matrices)
  invariantMatrices <- list(bmat, cmat, dmat, qmat,
                            rmat, xmat, pmat, umat)
  
  #### 4) create OpenMx models ####
  # create a list of models (one for each individual) for each latent class:
  model_list_k1 <- list()
  for(i in unique_ids){
    model_list_k1[[i]] <- mxModel(name=paste0('id_', i, "_k1"),
                                  amat_k1, bmat, cmat, dmat, qmat, rmat, xmat, pmat, umat,
                                  mxExpectationStateSpace('A1', 'B', 'C', 'D', 'Q', 'R', 'x0', 'P0', 'u'),
                                  mxFitFunctionML(),
                                  mxData(data[data[, id] == i, factors_ind], 'raw'))
  }
  
  model_list_k2 <- list()
  for(i in unique_ids){
    model_list_k2[[i]] <- mxModel(name=paste0('id_', i, "_k2"),
                                  amat_k2, bmat, cmat, dmat, qmat, rmat, xmat, pmat, umat,
                                  mxExpectationStateSpace('A2', 'B', 'C', 'D', 'Q', 'R', 'x0', 'P0', 'u'),
                                  mxFitFunctionML(),
                                  mxData(data[data[, id] == i, factors_ind], 'raw'))
  }
  
  if(n_clusters == 4){
    model_list_k3 <- list()
    for(i in unique_ids){
      model_list_k3[[i]] <- mxModel(name=paste0('id_', i, "_k3"),
                                    amat_k3, bmat, cmat, dmat, qmat, rmat, xmat, pmat, umat,
                                    mxExpectationStateSpace('A3', 'B', 'C', 'D', 'Q', 'R', 'x0', 'P0', 'u'),
                                    mxFitFunctionML(),
                                    mxData(data[data[, id] == i, factors_ind], 'raw'))
    }
    
    model_list_k4 <- list()
    for(i in unique_ids){
      model_list_k4[[i]] <- mxModel(name=paste0('id_', i, "_k4"),
                                    amat_k4, bmat, cmat, dmat, qmat, rmat, xmat, pmat, umat,
                                    mxExpectationStateSpace('A4', 'B', 'C', 'D', 'Q', 'R', 'x0', 'P0', 'u'),
                                    mxFitFunctionML(),
                                    mxData(data[data[, id] == i, factors_ind], 'raw'))
    }
  }
  
  modelnames_k1 <- paste0('id_', unique_ids, "_k1")
  modelnames_k2 <- paste0('id_', unique_ids, "_k2")
  if(n_clusters == 4){
    modelnames_k3 <- paste0('id_', unique_ids, "_k3")
    modelnames_k4 <- paste0('id_', unique_ids, "_k4")
  }
  
  # vector of model names (needed later)
  objectives_k1 <- paste(modelnames_k1, "objective", sep = ".")
  objectives_k2 <- paste(modelnames_k2, "objective", sep = ".")
  if(n_clusters == 4){
    objectives_k3 <- paste(modelnames_k3, "objective", sep = ".")
    objectives_k4 <- paste(modelnames_k4, "objective", sep = ".")
  }
  # vector of names of the objective functions (needed later)
  
  #### 5) mixture modeling ####
  # maximum number of iterations:
  #set.seed(8389493)#!!!! work on the whole replicability thing!!!!
  # provide seeds for the multiple starts (for replicability)
  seeds <- sample(1:100000000, nstarts)
  nonconvergences <- 0
  estimation_start <- Sys.time()
  # loop across the random starts
  for(i in 1:nstarts){
    # set seed:
    set.seed(seeds[i])
    # create random cluster assignment to start:
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
    
    # loop over iterations until convergence or max iterations are reached
    for(it in 1:maxit) {
      ## M-step: fitting SSM model and update class proportions
      
      # create a multi-person model for each latent class
      # the fit function is weighted by the posterior probabilities of each person-model
      # create the mxModel for cluster 1:
      weights = post[, 1]
      weighted_objectives <- paste(weights, "*", objectives_k1, collapse = " + ")
      model_k1 <- mxModel('model_k1', model_list_k1, mxAlgebraFromString(weighted_objectives, name = "weightedfit"), mxFitFunctionAlgebra("weightedfit"))
      
      # create the mxModel for cluster 2:
      weights = post[, 2]
      weighted_objectives <- paste(weights, "*", objectives_k2, collapse = " + ")
      model_k2 <- mxModel('model_k2', model_list_k2, mxAlgebraFromString(weighted_objectives, name = "weightedfit"), mxFitFunctionAlgebra("weightedfit"))
      
      # create models for clusters 3 and 4 (if applicable):
      if(n_clusters == 4){
        # create the mxModel for cluster 3:
        weights = post[, 3]
        weighted_objectives <- paste(weights, "*", objectives_k3, collapse = " + ")
        model_k3 <- mxModel('model_k3', model_list_k3, mxAlgebraFromString(weighted_objectives, name = "weightedfit"), mxFitFunctionAlgebra("weightedfit"))
        
        # create the mxModel for cluster 4:
        weights = post[, 4]
        weighted_objectives <- paste(weights, "*", objectives_k4, collapse = " + ")
        model_k4 <- mxModel('model_k4', model_list_k4, mxAlgebraFromString(weighted_objectives, name = "weightedfit"), mxFitFunctionAlgebra("weightedfit"))
      }
      
      # combine all models in one full model
      if(n_clusters == 2){
        fullmodel <- mxModel("fullmodel", model_k1, model_k2, mxFitFunctionMultigroup(c("model_k1", "model_k2")))
      }
      if(n_clusters == 4){
        fullmodel <- mxModel("fullmodel", model_k1, model_k2, model_k3, model_k4,
                             mxFitFunctionMultigroup(c("model_k1", "model_k2", "model_k3", "model_k4"))
        )
      }
      
      
      fullmodelr <- mxRun(fullmodel, silent = !verbose, suppressWarnings = TRUE)
      
      # compute class proportions:
      class_proportions <- colMeans(post)
      
      ## E-step
      # computing casewise log likelihoods in all models and total loglik
      fitfunctions_k1 <- paste0(modelnames_k1, ".fitfunction")
      fitfunctions_k2 <- paste0(modelnames_k2, ".fitfunction")
      if(n_clusters == 4){
        fitfunctions_k3 <- paste0(modelnames_k3, ".fitfunction")
        fitfunctions_k4 <- paste0(modelnames_k4, ".fitfunction")
      }
      # obtain the minus-2 log likelihoods from the OpenMx output
      m2_loglik_k1 <- purrr::map_dbl(fitfunctions_k1, mxEvalByName, model=fullmodelr)
      m2_loglik_k2 <- purrr::map_dbl(fitfunctions_k2, mxEvalByName, model=fullmodelr)
      if(n_clusters == 4){
        m2_loglik_k3 <- purrr::map_dbl(fitfunctions_k3, mxEvalByName, model=fullmodelr)
        m2_loglik_k4 <- purrr::map_dbl(fitfunctions_k4, mxEvalByName, model=fullmodelr)
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
      
      # compute total log likelihood
      total_loglik = sum(casewise.loglik*post)
      
      
      if(verbose){
        if(it != 1){
          print(paste0("Iteration: ", it, ". Log Likelihood: ", round(total_loglik, 4), ". Change: ", round(total_loglik - loglik0, 6), "."))
        } else {
          print(paste0("Iteration: ", it, ". Log Likelihood: ", round(total_loglik, 4), "."))
        }
        
      }
      # check convergence:
      if(it>1 && (total_loglik - loglik0)<1.0e-6){
        if(verbose){
          print("Convergence achieved.")
        }
        if(total_loglik < loglik0){
          fullmodelr <- fullmodelr0
          total_loglik <- loglik0
        }
        break
      }
      # update posteriors:
      post <- EStep(pi = class_proportions, ngroup = n, nclus = n_clusters, loglik = casewise.loglik)
      
      loglik0 = total_loglik
      fullmodelr0 = fullmodelr
      
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
      best_loglik <- total_loglik
      best_post <- post
      most_iterations <- it
      best_model <- fullmodelr
    } else {
      # otherwise, compare fits
      if(total_loglik > best_loglik){
        best_loglik <- total_loglik
        best_post <- post
        best_model <- fullmodelr
      }
      if(it > most_iterations){
        most_iterations <- it
      }
    }
    
    if(verbose){
      print(paste("Start", i, "out of", nstarts, "completed."))
    }
  }
  duration <- difftime(Sys.time(), estimation_start, unit = "s")
  
  #### find proxy maximum ###
  # !!! TO BE ADDED !!!
  
  #### 5) extract estimates ####
  ## TO DO: make this more flexible 
  # cluster 1:
  estimates <- coef(best_model)
  
  loglik <- best_loglik
  post <- best_post
  assignment <- round(post)
  class_proportions <- colMeans(post)
  
  clustering <- list("class_proportions" = class_proportions,
                     "posterior_prob" = post,
                     "assignment" = assignment)
  
  other <- list("loglik" = loglik,
                "nonconvergences" = nonconvergences,
                #"proxy_maximum" = proxy_maximum, #remove for package
                "duration" = duration,
                "most_iterations" = most_iterations)
  
  #### 5) build the output ####
  output <- list("data" = data,
                 "estimates" = estimates,
                 "clustering" = clustering,
                 "other" = other)
  
  return(output)
}