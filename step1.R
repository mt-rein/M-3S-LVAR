step1 <- function(data, measurementmodel, id,
                  invariances = c("loadings", "intercepts"),
                  partial_noninvariances = NULL){
  # data:
  #   a data frame with the indicator and ID variables
  # measurementmodel:
  #   a string describing the measurement model using the lavaan syntax
  # id:
  #   a character that indicates the id variable (the variable that indicates
  #   which observations belong to which person)
  
  # Are there measurement blocks in the MM?
  if(is.list(measurementmodel)){
    # If there are measurement blocks, how many do we have?
    M <- length(measurementmodel)
    
    MMoutput <- vector("list", length = M)
    for(m in 1:M){
      # estimate MM in each block
      MMoutput[[m]] <- lavaan::cfa(measurementmodel[[m]],
                                   data = data,
                                   estimator = "ML",
                                   group = id,
                                   group.equal = invariances,
                                   group.partial = partial_noninvariances,
                                   se = "none", test = "none",
                                   baseline = FALSE, h1 = FALSE)
    }
  } else {
    # if there are no measurement blocks, estimate full MM
    MMoutput <- lavaan::cfa(measurementmodel,
                            data = data,
                            estimator = "ML",
                            group = id,
                            group.equal = invariances,
                            group.partial = partial_noninvariances,
                            se = "none", test = "none",
                            baseline = FALSE, h1 = FALSE)
  }
  
  # assemble output
  output <- list("MMoutput" = MMoutput,
                 "data" = data,
                 "measurementmodel" = measurementmodel)
  return(output)
}