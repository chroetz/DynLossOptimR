predictDefault <- function(model, initialConditions, nPred) {
  d <- ncol(initialConditions)
  prediction <- array(NA_real_, dim = c(nPred, d, nrow(initialConditions)))
  initialConditions <- PolyPropR::normalize(initialConditions, model$normalization)
  currentState <- initialConditions
  for (i in seq_len(nPred)) {
    currentState <- PolyPropR::evaluateMonomialFeatures(currentState, model$nDeg) %*% model$coef
    prediction[i, , ] <- PolyPropR::denormalize(currentState, model$normalization)
  }
  return(prediction)
}
