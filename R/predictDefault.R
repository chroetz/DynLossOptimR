predictDefault <- function(model, initialConditions, nPred, softSquashScale=1e4) {
  d <- ncol(initialConditions)
  prediction <- array(NA_real_, dim = c(nPred, d, nrow(initialConditions)))
  initialConditions <- PolyPropR::normalize(initialConditions, model$normalization)
  currentState <- initialConditions
  for (i in seq_len(nPred)) {
    for (j in seq_len(model$intermediate)) {
        currentState <- PolyPropR::evaluateMonomialFeatures(currentState, model$nDeg) %*% model$coef
        currentState <- softSquashScale * tanh(currentState / softSquashScale)
    }
    prediction[i, , ] <- t(PolyPropR::denormalize(currentState, model$normalization))
  }
  return(prediction)
}
