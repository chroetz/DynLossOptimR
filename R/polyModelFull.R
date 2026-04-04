#' @export
fitFullLoss <- function(
  xTrain,
  nDeg,
  weightSchedule,
  normalizationType = "none",
  normalizationScale = 1.0,
  targetUpdateFactor = 0
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
      data = c(model = "poly_model_full",
               list(
                 obs = target,
                 weights_obs = ws$weightsObs,
                 weights_ana = ws$weightsAna,
                 weights_pen = ws$weightsPen,
                 weights_coef = ws$weightsCoef,
                 deg = nDeg
               )),
      parameters = list(coef = coef, ana = analysis),
      DLL = "DynLossOptimR_TMBExports",
      silent = TRUE
    )

    opt <- optimizer(obj$par, obj$fn, obj$gr, reltol = ws$reltol)

    analysisNew <- opt$par[names(opt$par) == "ana"]
    coefNew <- opt$par[names(opt$par) == "coef"]
    dim(analysisNew) <- dim(analysis)
    dim(coefNew) <- dim(coef)

    analysis <- analysisNew
    coef <- coefNew
    target <- (1 - targetUpdateFactor) * target + targetUpdateFactor * analysis
  }

  return(
    list(
      coef = coef,
      nDeg = nDeg,
      analysis = PolyPropR::denormalize(analysis, normalization),
      normalization = normalization
    )
  )

}

predictFullLoss <- function(model, initialConditions, nPred) {
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
