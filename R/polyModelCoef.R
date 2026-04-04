#' @export
fitCoefLoss <- function(
  xTrain,
  nDeg,
  weightSchedule,
  normalizationType = "none",
  normalizationScale = 1.0
) {

  n <- nrow(xTrain)
  d <- ncol(xTrain)

  normalization <- PolyPropR::getNormalization(xTrain, normalizationType, normalizationScale)
  xTrain <- PolyPropR::normalize(xTrain, normalization)

  optimizer <- function(par, fn, gr, reltol = 1e-6) {
    stats::optim(par, fn, gr,
      method = "BFGS",
      control = list(
        maxit = 1e5,
        reltol = reltol
      )
    )
  }

  # Initial linear fit
  input1 <- xTrain[-n, , drop = FALSE]
  output1 <- xTrain[-1, , drop = FALSE]
  features1 <- PolyPropR::evaluateMonomialFeatures(input1, nDeg)
  coef <- PolyPropR::fitLinear(features1, output1)

  analysis <- xTrain
  target <- xTrain

  for (k in seq_along(weightSchedule)) {
    ws <- weightSchedule[[k]]

    obj <- TMB::MakeADFun(
      data = c(model = "poly_model_coef",
               list(
                 obs = target,
                 weights_obs = ws$weightsObs,
                 weights_coef = ws$weightsCoef,
                 deg = nDeg
               )),
      parameters = list(coef = coef),
      DLL = "DynLossOptimR_TMBExports",
      silent = TRUE
    )

    opt <- optimizer(obj$par, obj$fn, obj$gr, reltol = ws$reltol)

    coefNew <- opt$par[names(opt$par) == "coef"]
    dim(coefNew) <- dim(coef)

    coef <- coefNew
  }

  return(
    list(
      coef = coef,
      nDeg = nDeg,
      analysis = PolyPropR::denormalize(xTrain, normalization),
      normalization = normalization
    )
  )

}

predictCoefLoss <- function(model, initialConditions, nPred) {
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
