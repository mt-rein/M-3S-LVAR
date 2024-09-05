#### This script defines the simulation function do_sim() ####

do_sim <- function(pos, cond, outputfile, verbose = FALSE){
  # pos = position in the condition grid
  # cond = the condition grid
  # outputfile = file name for the output CSV file
  # verbose = if TRUE, prints a message after the iteration is finished
  
  #### for testing:
  # pos = 1
  # replication <- 1
  # iteration <- 1
  # # get condition levels and set seed:
  # n <- 48
  # obs <- 50
  # n_k = 2
  # k_size =  "balanced" |> as.character()
  # rho_gen =  "high" |> as.character()
  # cluster_separation =  "high" |> as.character()
  # innovars =  "invariant" |> as.character()

  
  replication <- cond$replication[pos]
  iteration <- cond$iteration[pos]
  # get condition levels and set seed:
  n <- cond$n[pos]
  obs <- cond$obs[pos]
  n_k = cond$n_k[pos]
  k_size =  cond$k_size[pos] |> as.character()
  rho_gen =  cond$rho_gen[pos] |> as.character()
  cluster_separation =  cond$cluster_separation[pos] |> as.character()
  innovars =  cond$innovars[pos] |> as.character()
  seed_cond <- cond$seed[pos]
  set.seed(seed_cond)
  
  # change some global options in OpenMx:
  mxOption(key = "Calculate Hessian", value = "No")
  mxOption(key = "Standard Errors", value = "No") 
  
  #### set data generation parameters ####
  ## regression parameters:
  # cluster 1
  if(cluster_separation == "low"){
    phi11_k1_pop <- .3
    phi22_k1_pop <- .6
    phi12_k1_pop <- .3
    phi21_k1_pop <- .3
  }
  if(cluster_separation == "high"){
    phi11_k1_pop <- 0
    phi22_k1_pop <- .6
    phi12_k1_pop <- .3
    phi21_k1_pop <- .3
  }
  
  # cluster 2
  if(cluster_separation == "low"){
    phi11_k2_pop <- .6
    phi22_k2_pop <- .6
    phi12_k2_pop <- 0
    phi21_k2_pop <- .3
  }
  if(cluster_separation == "high"){
    phi11_k2_pop <- .6
    phi22_k2_pop <- .6
    phi12_k2_pop <- -.3
    phi21_k2_pop <- .3
  }
  
  # clusters 3 and 4
  if(n_k  == 4){
    ## in case of 4 clusters, set regression parameters
    # cluster 3:
    if(cluster_separation == "low"){
      phi11_k3_pop <- .6
      phi22_k3_pop <- .3
      phi12_k3_pop <- .3
      phi21_k3_pop <- .3
    }
    if(cluster_separation == "high"){
      phi11_k3_pop <- .6
      phi22_k3_pop <- 0
      phi12_k3_pop <- .3
      phi21_k3_pop <- .3
    }
    
    # cluster 4:
    if(cluster_separation == "low"){
      phi11_k4_pop <- .6
      phi22_k4_pop <- .6
      phi12_k4_pop <- .3
      phi21_k4_pop <- 0
    }
    if(cluster_separation == "high"){
      phi11_k4_pop <- .6
      phi22_k4_pop <- .6
      phi12_k4_pop <- .3
      phi21_k4_pop <- -.3
    }
    
  } else {
    ## in case of 2 clusters, set to NA
    # cluster 3:
    phi11_k3_pop <- NA
    phi22_k3_pop <- NA
    phi12_k3_pop <- NA
    phi21_k3_pop <- NA
    
    # cluster 4:
    phi11_k4_pop <- NA
    phi22_k4_pop <- NA
    phi12_k4_pop <- NA
    phi21_k4_pop <- NA
  }
  
  # combine into matrices:
  phimat_k1 <- matrix(c(phi11_k1_pop, phi21_k1_pop,
                        phi12_k1_pop, phi22_k1_pop),
                      ncol = 2)
  phimat_k2 <- matrix(c(phi11_k2_pop, phi21_k2_pop,
                        phi12_k2_pop, phi22_k2_pop),
                      ncol = 2)
  if(n_k == 4){
    phimat_k3 <- matrix(c(phi11_k3_pop, phi21_k3_pop,
                          phi12_k3_pop, phi22_k3_pop),
                        ncol = 2)
    phimat_k4 <- matrix(c(phi11_k4_pop, phi21_k4_pop,
                          phi12_k4_pop, phi22_k4_pop),
                        ncol = 2)
  }
  
  
  ## innovation variances
  # fixed effect innovation covariance matrix is invariant across clusters
  zeta1_pop <- 1.5
  zeta2_pop <- 1.5
  zeta12_pop <- .5
  
  ## means
  # cluster mean is invariant across clusters
  grandmeans <- c(5, 5)
  
  #### generate factor scores ####
  # create empty data frame
  eta <- tibble(id = numeric(),
                obs = numeric(),
                eta1 = numeric(),
                eta2 = numeric(),
                k_true = numeric())
  
  ## create cluster assignment vector
  # for two clusters:
  if(n_k == 2){
    # for balanced distribution:
    if(k_size == "balanced"){
      clusterassignment_true <- c(rep(1, n*0.5), rep(2, n*0.5))
    }
    # for unbalanced distribution
    if(k_size == "unbalanced"){
      clusterassignment_true <- c(rep(1, n*0.75), rep(2, n*0.25))
    }
  }
  if(n_k == 4){
    # for balanced distribution:
    if(k_size == "balanced"){
      clusterassignment_true <- c(rep(1, n*0.25), rep(2, n*0.25), rep(3, n*0.25), rep(4, n*0.25))
    }
    # for unbalanced distribution
    if(k_size == "unbalanced"){
      clusterassignment_true <- c(rep(1, n*0.75), rep(2, n*0.25/3), rep(3, n*0.25/3), rep(4, n*0.25/3))
    }
  }
  
  for(i in 1:n){
    # create temporary dataframe:
    k <- clusterassignment_true[i]
    
    ## create (partly person-specific) data-generating parameter values from population values
    # get correct phi matrix:
    if(k == 1){
      phimat <- phimat_k1
    }
    if(k == 2){
      phimat <- phimat_k2
    }
    if(k == 3){
      phimat <- phimat_k3
    }
    if(k == 4){
      phimat <- phimat_k4
    }
    
    # set innovation variance matrix:
    if(innovars == "invariant"){
      # if innovation (co)variance matrix is invariant across persons, set to population values
      zeta1_i <- zeta1_pop
      zeta2_i <- zeta2_pop
      zeta12_i <- zeta12_pop
    }
    if(innovars == "random"){
      # otherwise, create person specific innovation (co)variance matrix:
      zeta1_i <- rtruncnorm(1, a = zeta1_pop - .5, b = zeta1_pop + .5,
                            mean = zeta1_pop, sd = .5)
      zeta2_i <- rtruncnorm(1, a = zeta2_pop - .5, b = zeta2_pop + .5,
                            mean = zeta2_pop, sd = .5)
      zeta12_i <- rtruncnorm(1, a = zeta12_pop - .5, b = zeta12_pop + .5,
                             mean = zeta12_pop, sd = .5)
    }
    zetamat <- matrix(c(zeta1_i, zeta12_i, zeta12_i, zeta2_i), ncol = 2)
    
    mu_i <- c(rnorm(1, mean = grandmeans[1], sd = 2),                           # person specific mean on variable 1
              rnorm(1, mean = grandmeans[2], sd = 2))                           # person specific mean on variable 2
    
    eta_i <- sim_VAR(factors = 2, obs = obs,
                     phi = phimat, zeta = zetamat,
                     mu = mu_i,
                     burn_in = 10)
    
    # add id and true cluster variable:
    eta_i$id <- i
    eta_i$k_true <- k
    
    # merge with full data
    eta <- dplyr::full_join(eta, eta_i, by = join_by(id, obs, eta1, eta2, k_true))
  }
  
  #### generate indicators ####
  # lambda matrix (loadings), all 1:
  loadings <- rep(1, 4)
  lambda <- matrix(c(loadings, rep(0, 4), rep(0, 4), loadings),
                   nrow = 8, ncol = 2)
  
  # get "true" psi (factor variance)
  zetamat <- matrix(c(zeta1_pop, zeta12_pop, zeta12_pop, zeta2_pop), ncol = 2)
  psi_k1 <- solve(diag(2*2) - kronecker(phimat_k1, phimat_k1)) %*% c(zetamat) |> matrix(nrow = 2)
  psi_k2 <- solve(diag(2*2) - kronecker(phimat_k2, phimat_k2)) %*% c(zetamat) |> matrix(nrow = 2)
  if(n_k == 4){
    psi_k3 <- solve(diag(2*2) - kronecker(phimat_k3, phimat_k3)) %*% c(zetamat) |> matrix(nrow = 2)
    psi_k4 <- solve(diag(2*2) - kronecker(phimat_k4, phimat_k4)) %*% c(zetamat) |> matrix(nrow = 2)
  }
  
  n_k1 <- sum(clusterassignment_true == 1)
  n_k2 <- sum(clusterassignment_true == 2)
  
  if(n_k == 2){
    psi_pop <- (psi_k1*(n_k1-1) + psi_k2*n_k2)/(n_k1 + n_k2 - 2)
  }
  if(n_k == 4){
    n_k3 <- sum(clusterassignment_true == 3)
    n_k4 <- sum(clusterassignment_true == 4)
    psi_pop <- (psi_k1 + psi_k2 + psi_k3 + psi_k4)/(n_k1 + n_k2 + n_k3 + n_k4 - 4)
  }
  
  
  # write an objective function that gives the difference between the computed rho
  # (based on the error variances) and desired rho
  rho_difference <- function(errorvars, psi, lambda, target_rho) {
    theta <- matrix(0, nrow = 8, ncol = 8)                                      # error variance matrix with 0 on off-diagonal
    diag(theta) <- errorvars                                                    # add error variances to the diagonal of theta
    computed_rho <- diag(psi %*% t(lambda) %*% solve(lambda %*% psi %*% t(lambda) + theta) %*% lambda)
    diff_rho <- sum((computed_rho - target_rho)^2)                              # squared difference between computed and desired rho as the objective
    return(diff_rho)
  }
  
  # initial values for error variances
  initial_errorvars <- 0.5
  
  # target values for rho
  if(rho_gen == "low"){
    target_rho <- c(0.6, 0.6)
  }
  if(rho_gen == "high"){
    target_rho <- c(0.9, 0.9)
  }
  
  
  # run optimization
  result <- optim(
    par = initial_errorvars,
    fn = rho_difference,
    psi = psi_pop,  # Provide psi, lambda, and target_rho as additional arguments
    lambda = lambda,
    target_rho = target_rho,
    method = "L-BFGS-B",
    lower = rep(0.0001, 8)
  )
  
  # retrieve the optimized theta
  optimized_errorvars <- result$par
  
  theta <- matrix(0, nrow = 8, ncol = 8)                                        # create error covariance matrix
  diag(theta) <- optimized_errorvars                                            # place error variances on theta diagonal
  
  epsilon <- mvrnorm(nrow(eta), mu = rep(0, 8), Sigma = theta, empirical=T)     # generate errors
  # generate indicator scores
  # (intercepts left out because we set them to 0):
  data <- as.matrix(eta[, c("eta1", "eta2")]) %*% t(lambda) + epsilon %>%
    as.data.frame()
  colnames(data) <- paste0("v", 1:8)
  data$id <- eta$id
  data$obs <- eta$obs
  data$k_true <- eta$k_true
  
  
  #### Step 1 ####
  model_step1 <- list("f1 =~ v1 + v2 + v3 + v4",
                      "f2 =~ v5 + v6 + v7 + v8")
  
  output_step1 <- run_step1(data = data,
                            measurementmodel = model_step1,
                            id = "id",
                            invariances = c("loadings", "intercepts"),
                            partial_noninvariances = NULL)
  # extract error/warning messages (if applicable):
  step1_warning <- ifelse(is_empty(output_step1$warnings),
                          FALSE, TRUE)
  step1_warning_text <- ifelse(is_empty(output_step1$warnings),
                               "",
                               paste(c(output_step1$warnings),
                                     collapse = "; ")
  )
  step1_error <- ifelse(is_empty(output_step1$result$error),
                        FALSE, TRUE)
  step1_error_text <- ifelse(is_empty(output_step1$result$error),
                             "",
                             paste(c(output_step1$result$error),
                                   collapse = "; "))
  
  #### Step 2 ####
  if(!step1_error){                                                             # only proceed if there is no error in step 1
    output_step2 <- run_step2(step1output = output_step1$result$result)
    # extract error/warning messages (if applicable):
    step2_warning <- ifelse(is_empty(output_step2$warnings),
                            FALSE, TRUE)
    step2_warning_text <- ifelse(is_empty(output_step2$warnings),
                                 "",
                                 paste(c(output_step2$warnings),
                                       collapse = "; ")
    )
    step2_error <- ifelse(is_empty(output_step2$result$error),
                          FALSE, TRUE)
    step2_error_text <- ifelse(is_empty(output_step2$result$error),
                               "",
                               paste(c(output_step2$result$error),
                                     collapse = "; ")
    )
    # lambdas_est <- lavInspect(output_step1$result$result$fit_step1, "est")$lambda
    # lambda2 <- lambdas_est[2,1]
    # lambda3 <- lambdas_est[3,1]
    # lambda4 <- lambdas_est[4,1]
    # lambda5 <- lambdas_est[5,2]
    # lambda6 <- lambdas_est[6,1]
    # lambda7 <- lambdas_est[7,1]
    # lambda8 <- lambdas_est[8,1]
    # 
    # thetas_est <- lavInspect(output_step1$result$result$fit_step1, "est")$theta
    # theta1 <- thetas_est[1,1]
    # theta2 <- thetas_est[2,2]
    # theta3 <- thetas_est[3,3]
    # theta4 <- thetas_est[4,4]
    # theta5 <- thetas_est[5,5]
    # theta6 <- thetas_est[6,6]
    # theta7 <- thetas_est[7,7]
    # theta8 <- thetas_est[8,8] 
  } else {
    step2_warning <- FALSE
    step2_warning_text <- "step1 not successful"
    step2_error <- FALSE
    step2_error_text <- "step1 not successful"
    
    # lambda1 <- NA
    # lambda2 <- NA
    # lambda3 <- NA
    # lambda4 <- NA
    # lambda5 <- NA
    # lambda6 <- NA
    # lambda7 <- NA
    # lambda8 <- NA
    # 
    # theta1 <- NA
    # theta2 <- NA
    # theta3 <- NA
    # theta4 <- NA
    # theta5 <- NA
    # theta6 <- NA
    # theta7 <- NA
    # theta8 <- NA 
  }
  
  #### Step 3 ####
  if(!step1_error & !step2_error){                                              # only proceed if there is no error in step 1 as well as step 2
    output_step3 <- run_step3(step2output = output_step2$result$result,
                              n_clusters = n_k,
                              n_starts = 15, n_best_starts = 5,
                              maxit = 100,
                              true_clusters = clusterassignment_true, 
                              verbose = FALSE)
    
    # extract error/warning messages (if applicable):
    step3_warning <- ifelse(is_empty(output_step3$warnings),
                            FALSE, TRUE)
    step3_warning_text <- ifelse(is_empty(output_step3$warnings),
                                 "",
                                 paste(c(output_step3$warnings),
                                       collapse = "; ")
    )
    step3_error <- ifelse(is_empty(output_step3$result$error),
                          FALSE, TRUE)
    step3_error_text <- ifelse(is_empty(output_step3$result$error),
                               "",
                               paste(c(output_step3$result$error),
                                     collapse = "; ")
    )
  } else {
    step3_warning <- FALSE
    step3_warning_text <- "step1 or step2 not successful"
    step3_error <- FALSE
    step3_error_text <- "step1 or step2 not successful"
  }
  
  if(!step1_error & !step2_error & !step3_error){
    final_output <- output_step3$result$result
    duration = final_output$other$duration |> as.numeric()
    nonconvergences = final_output$other$nonconvergences
    
    ## adjust potential label switching:
    combinations <- RcppAlgos::permuteGeneral(paste0("cluster", 1:n_k)) |> as.data.frame()
    combinations$diagsum <- 0
    
    for(i in 1:nrow(combinations)){
      # relabel the clusters:
      colnames(final_output$clustering$modal_assignment) <- combinations[i, 1:n_k]
      clusterassignment_estimated = character(n)
      for(j in 1:n){
        index <- which(final_output$clustering$modal_assignment[j, ] == 1)
        clusterassignment_estimated[j] <- colnames(final_output$clustering$modal_assignment)[index]
      }
      
      # creates a cross table of estimated and true cluster assignments:
      crosstable <- table(clusterassignment_estimated, clusterassignment_true)
      # compute the sum of the diagonal of the cross table and save it
      combinations$diagsum[i] <- sum(diag(crosstable))                          
    }
    # choose the label with the largest sum of the diagonal:
    newlabels <- combinations[which.max(combinations$diagsum), 1:n_k]
    
    # swap the labels accordingly in the output of step 3:
    names(final_output$estimates) <- colnames(final_output$clustering$posterior_prob) <- names(final_output$clustering$class_proportions) <- colnames(final_output$clustering$modal_assignment) <- newlabels
    
    ## extract estimates
    # estimates in cluster 1:
    phi11_k1 <- final_output$estimates$cluster1["phi11"] |> as.numeric()
    phi12_k1 <- final_output$estimates$cluster1["phi12"] |> as.numeric()
    phi21_k1 <- final_output$estimates$cluster1["phi21"] |> as.numeric()
    phi22_k1 <- final_output$estimates$cluster1["phi22"] |> as.numeric()
    
    zeta1_k1 <- final_output$estimates$cluster1["zeta1"] |> as.numeric()
    zeta2_k1 <- final_output$estimates$cluster1["zeta2"] |> as.numeric()
    zeta12_k1 <- final_output$estimates$cluster1["zeta12"] |> as.numeric()
    
    # estimates in cluster 2:
    phi11_k2 <- final_output$estimates$cluster2["phi11"] |> as.numeric()
    phi12_k2 <- final_output$estimates$cluster2["phi12"] |> as.numeric()
    phi21_k2 <- final_output$estimates$cluster2["phi21"] |> as.numeric()
    phi22_k2 <- final_output$estimates$cluster2["phi22"] |> as.numeric()
    
    zeta1_k2 <- final_output$estimates$cluster2["zeta1"] |> as.numeric()
    zeta2_k2 <- final_output$estimates$cluster2["zeta2"] |> as.numeric()
    zeta12_k2 <- final_output$estimates$cluster2["zeta12"] |> as.numeric()
    
    # estimates in cluster 3 and 4 (if applicable)
    if(n_k == 4){
      # estimates in cluster 3:
      phi11_k3 <- final_output$estimates$cluster3["phi11"] |> as.numeric()
      phi12_k3 <- final_output$estimates$cluster3["phi12"] |> as.numeric()
      phi21_k3 <- final_output$estimates$cluster3["phi21"] |> as.numeric()
      phi22_k3 <- final_output$estimates$cluster3["phi22"] |> as.numeric()
      
      zeta1_k3 <- final_output$estimates$cluster3["zeta1"] |> as.numeric()
      zeta2_k3 <- final_output$estimates$cluster3["zeta2"] |> as.numeric()
      zeta12_k3 <- final_output$estimates$cluster3["zeta12"] |> as.numeric()
      
      # estimates in cluster 4:
      phi11_k4 <- final_output$estimates$cluster4["phi11"] |> as.numeric()
      phi12_k4 <- final_output$estimates$cluster4["phi12"] |> as.numeric()
      phi21_k4 <- final_output$estimates$cluster4["phi21"] |> as.numeric()
      phi22_k4 <- final_output$estimates$cluster4["phi22"] |> as.numeric()
      
      zeta1_k4 <- final_output$estimates$cluster4["zeta1"] |> as.numeric()
      zeta2_k4 <- final_output$estimates$cluster4["zeta2"] |> as.numeric()
      zeta12_k4 <- final_output$estimates$cluster4["zeta12"] |> as.numeric()
      
      
    } else {
      # set everything to NA if only 2 clusters
      phi11_k3 <- NA
      phi12_k3 <- NA
      phi21_k3 <- NA
      phi22_k3 <- NA
      
      zeta1_k3 <- NA
      zeta2_k3 <- NA
      zeta12_k3 <- NA
      
      phi11_k4 <- NA
      phi12_k4 <- NA
      phi21_k4 <- NA
      phi22_k4 <- NA
      
      zeta1_k4 <- NA
      zeta2_k4 <- NA
      zeta12_k4 <- NA
    }
    
    # ARI:
    clusterassignment_estimated <- apply(final_output$clustering$modal_assignment, 1, 
                                         function(row) {
                                           colnames(final_output$clustering$modal_assignment)[which(row == 1)]
                                         })
    ARI <- mcclust::arandi(clusterassignment_true, clusterassignment_estimated, adjust = TRUE)
    
    # local maximum:
    local_max <- abs(final_output$other$proxy_maximum - final_output$other$loglik) > .001
  } else {
    duration = NA
    nonconvergences = NA
    phi11_k1 <- NA
    phi12_k1 <- NA
    phi21_k1 <- NA
    phi22_k1 <- NA
    
    zeta1_k1 <- NA
    zeta2_k1 <- NA
    zeta12_k1 <- NA
    
    phi11_k2 <- NA
    phi12_k2 <- NA
    phi21_k2 <- NA
    phi22_k2 <- NA
    
    zeta1_k2 <- NA
    zeta2_k2 <- NA
    zeta12_k2 <- NA
    
    phi11_k3 <- NA
    phi12_k3 <- NA
    phi21_k3 <- NA
    phi22_k3 <- NA
    
    zeta1_k3 <- NA
    zeta2_k3 <- NA
    zeta12_k3 <- NA
    
    phi11_k4 <- NA
    phi12_k4 <- NA
    phi21_k4 <- NA
    phi22_k4 <- NA
    
    zeta1_k4 <- NA
    zeta2_k4 <- NA
    zeta12_k4 <- NA
    
    ARI <- NA
    local_max <- NA
  }
  
  # reset the OpenMx parameters:
  mxOption(key = "Calculate Hessian", reset = TRUE) 
  mxOption(key = "Standard Errors", reset = TRUE)
  mxOption(key = "Major Iterations", reset = TRUE)
  
  output <- c("iteration" = iteration, "replication" = replication,
              "n" = n, "obs" = obs, "n_k" = n_k, "k_size" = k_size, "rho_gen" = rho_gen, 
              "cluster_separation" = cluster_separation, "innovars" = innovars,
              "duration" = duration, "nonconvergences" = nonconvergences, "ARI" = ARI, "local_max" = local_max,
              "phi11_k1_pop" = phi11_k1_pop, "phi12_k1_pop" = phi12_k1_pop, "phi21_k1_pop" = phi21_k1_pop, "phi22_k1_pop" = phi22_k1_pop,
              "phi11_k2_pop" = phi11_k2_pop, "phi12_k2_pop" = phi12_k2_pop, "phi21_k2_pop" = phi21_k2_pop, "phi22_k2_pop" = phi22_k2_pop,
              "phi11_k3_pop" = phi11_k3_pop, "phi12_k3_pop" = phi12_k3_pop, "phi21_k3_pop" = phi21_k3_pop, "phi22_k3_pop" = phi22_k3_pop,
              "phi11_k4_pop" = phi11_k4_pop, "phi12_k4_pop" = phi12_k4_pop, "phi21_k4_pop" = phi21_k4_pop, "phi22_k4_pop" = phi22_k4_pop,
              "zeta1_pop" = zeta1_pop, "zeta2_pop" = zeta2_pop, "zeta12_pop" = zeta12_pop,
              "lambda2" = lambda2, "lambda3" = lambda3, "lambda4" = lambda4, 
              "lambda5" = lambda5, "lambda6" = lambda6, "lambda7" = lambda7, "lambda8" = lambda8,
              "theta_pop" = optimized_errorvars,
              "theta1" = theta1, "theta2" = theta2, "theta3" = theta3, "theta4" = theta4, 
              "theta5" = theta5, "theta6" = theta6, "theta7" = theta7, "theta8" = theta8, 
              "phi11_k1" = phi11_k1, "phi12_k1" = phi12_k1, "phi21_k1" = phi21_k1, "phi22_k1" = phi22_k1,
              "phi11_k2" = phi11_k2, "phi12_k2" = phi12_k2, "phi21_k2" = phi21_k2, "phi22_k2" = phi22_k2,
              "phi11_k3" = phi11_k3, "phi12_k3" = phi12_k3, "phi21_k3" = phi21_k3, "phi22_k3" = phi22_k3,
              "phi11_k4" = phi11_k4, "phi12_k4" = phi12_k4, "phi21_k4" = phi21_k4, "phi22_k4" = phi22_k4,
              "zeta1_k1" = zeta1_k1, "zeta2_k1" = zeta2_k1, "zeta12_k1" = zeta12_k1,
              "zeta1_k2" = zeta1_k2, "zeta2_k2" = zeta2_k2, "zeta12_k2" = zeta12_k2,
              "zeta1_k3" = zeta1_k3, "zeta2_k3" = zeta2_k3, "zeta12_k3" = zeta12_k3,
              "zeta1_k4" = zeta1_k4, "zeta2_k4" = zeta2_k4, "zeta12_k4" = zeta12_k4,
              "step1_warning" = step1_warning, "step2_warning" = step2_warning, "step3_warning" = step3_warning,
              "step1_error" = step1_error, "step2_error" = step2_error, "step3_error" = step3_error,
              "seed" = seed_cond, "pos" = pos,
              "step1_warning_text" = step1_warning_text, "step2_warning_text" = step2_warning_text, "step3_warning_text" = step3_warning_text,
              "step1_error_text" = step1_error_text, "step2_error_text" = step2_error_text, "step3_error_text" = step3_error_text)
  
  for(i in 69:74){
    output[i] <- str_squish(output[i])                                          # removes all whitespace and linebreaks from the error and warning strings
    output[i] <- gsub(",", "", output[i])                                       # removes all commata from error and warning strings (to prevent messing up the CSV file)
  }
  
  # check if file exists
  if(!file.exists(outputfile)){
    # if file does not yet exist
    write.table(t(output), file = outputfile, append = FALSE, quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)
  } else {
    # lock the file to prevent multiple processes accessing it simultaneously
    lock <- flock::lock(outputfile)
    write.table(t(output), file = outputfile, append = TRUE, quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)
    # unlock the file
    flock::unlock(lock)
  }
  
  if(verbose == TRUE){
    print(paste("Simulation", pos, "completed at", Sys.time()))                 # prints a message when a replication is done (as a sign that R did not crash)
  }
  
  return(output)
}