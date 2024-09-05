step2 <- function(step1output){
  # step1output:
  #   the object that was generated using the step1() function
  
  ## Preparations
  fit_step1 <- step1output$MMoutput
  data <- step1output$data
  measurementmodel <- step1output$measurementmodel
  indicators <- lavaan::lavNames(lavaan::lavaanify(measurementmodel), "ov")
  factors <- lavaan::lavNames(lavaan::lavaanify(measurementmodel), "lv")
  if(is.list(fit_step1)){
    id <- lavInspect(fit_step1[[1]], "group")
    M <- length(measurementmodel)
  } else {
    id <- lavInspect(fit_step1, "group")
  }
  if(!is.character(data[, id])){
    #warning("Note: The ID variable has been transformed to a character.")
    data[, id] <- as.character(data[, id])
  }
  n_persons <- length(unique(data[, id]))
  
  ## compute factor scores
  # Are there measurement blocks?
  if(is.list(measurementmodel)){
    for(m in 1:M){
      # if there are measurement blocks, compute factor scores in each block
      temp <- lavPredict(fit_step1[[m]], assemble = TRUE, append.data = TRUE)
      # and append to original data
      data <- dplyr::full_join(data, temp, by = c(id, indicators[indicators %in% colnames(temp)]))
      }
  } else {
    # if there are no measurement blocks, compute factor scores for all latent 
    # variables simultaneously
      temp <- lavPredict(fit_step1, assemble = TRUE, append.data = TRUE)
      data <- dplyr::full_join(data, temp, by = c(id, indicators))
    }
  
  ## compute lambda_star and theta_star
  # create empty matrices to store values.
  # One row per individual, one column per factor
  lambda_star <- theta_star <- matrix(NA, nrow = n_persons, ncol = length(factors))
  
  # Are there measurement blocks?
  if(is.list(measurementmodel)){
    # if yes, create lists to store MM parameter values per block
    # (lists with M elements, where each element is a list with n_person elements)
    psi_block <- lambda_block <- theta_block <- vector(mode = "list", length = M)
    for (m in 1:M){
      EST_block         <- lavaan::lavInspect(fit_step1[[m]], "est")
      psi_block[[m]] <- lapply(X = EST_block, "[[", "psi")
      lambda_block[[m]] <- lapply(X = EST_block, "[[", "lambda")
      theta_block[[m]]  <- lapply(X = EST_block, "[[", "theta")
    }
    
    # create lists to store MM parameter values per person
    # (i.e., combine the block-wise MM parameters per person)
    psi_person <- lambda_person <- theta_person <-  vector(mode = "list", length = n_persons)
    for(i in 1:n_persons){
      # put MM matrices of each person in the respective list
      for (m in 1:M) {
        psi_person[[i]][[m]] <- psi_block[[m]][[i]]
        lambda_person[[i]][[m]] <- lambda_block[[m]][[i]]
        theta_person[[i]][[m]]  <- theta_block[[m]][[i]]
      }
      # combine the matrices of the different measurement blocks into a single matrix
      psi_person[[i]] <- lavaan::lav_matrix_bdiag(psi_person[[i]])
      lambda_person[[i]] <- lavaan::lav_matrix_bdiag(lambda_person[[i]])
      theta_person[[i]]  <- lavaan::lav_matrix_bdiag(theta_person[[i]])
      
      # name the matrices' rows and columns
      rownames(psi_person[[i]]) <- colnames(psi_person[[i]]) <- factors
      rownames(lambda_person[[i]]) <- indicators
      colnames(lambda_person[[i]]) <- factors
      rownames(theta_person[[i]]) <- colnames(theta_person[[i]]) <- indicators
    }
  } else {
    # if there are no measurement blocks, simply extract the MM parameters per person
    psi_person <- lambda_person <- theta_person <-  vector(mode = "list", length = n_persons)
    for(i in 1:n_persons){
      psi_person[[i]] <- lavaan::lavInspect(fit_step1, "est")[[i]]$psi
      lambda_person[[i]] <- lavaan::lavInspect(fit_step1, "est")[[i]]$lambda
      theta_person[[i]] <- lavaan::lavInspect(fit_step1, "est")[[i]]$theta
    }
  }
  
  # compute lambda_star and theta_star per person
  for (i in 1:n_persons) {
    sigma_i <- lambda_person[[i]] %*% psi_person[[i]] %*% t(lambda_person[[i]]) + theta_person[[i]]
    lambda_star[i, ] <- diag(psi_person[[i]] %*% t(lambda_person[[i]]) %*% solve(sigma_i) %*% lambda_person[[i]])
    theta_star[i, ] <- lambda_star[i, ]*(1-lambda_star[i, ])*diag(psi_person[[i]])
  }
  
  # assemble output
  output <- list("data" = data,
                 "lambda_star" = lambda_star,
                 "theta_star" = theta_star,
                 "other" = list("factors" = factors,
                                "indicators" =  indicators,
                                "id" = id))
  return(output)
  }