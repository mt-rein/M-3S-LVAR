# only for testing!
# load("examples/outputstep1.RData")
# step1output <- output

step2 <- function(step1output){
  # step1output:
  #   the object that was generated using the step1() function
  
  ## Preparations
  fit_step1 <- step1output$fit_step1
  data <- step1output$data
  indicators <- lavaan::lavNames(fit_step1, "ov")
  factors <- lavaan::lavNames(fit_step1, "lv")
  
  ## compute factor scores:
  factorscores <- lavaan::lavPredict(fit_step1,
                                     assemble = TRUE,
                                     append.data = TRUE) |> 
    as.data.frame()
  # append factor scores to original data:
  data <- dplyr::full_join(data, factorscores, by = indicators)
  
  ## compute rho and kappa
  # These values will be used in step 3 to account for measurement error
  # rho indicates the model-based reliability:
  psi <- lavaan::lavInspect(fit_step1, what = "est")$psi
  lambda <- lavaan::lavInspect(fit_step1, what = "est")$lambda
  sigma <- fitted(fit_step1)$cov
  rho <- diag(psi %*% t(lambda) %*% solve(sigma) %*% lambda)
  
  # kappa indicates the uncertainty of the factor score estimation:
  kappa <- rho*(1-rho)*diag(psi)
  
  # assemble output
  output <- list("data" = data,
                 "rho" = rho,
                 "kappa" = kappa,
                 "fit_step1" = fit_step1)
  return(output)
}