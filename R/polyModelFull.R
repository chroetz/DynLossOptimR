#' @export
fitFullLoss <- function(
  xTrain,
  nDeg,
  weightSchedule,
  normalizationType = "none",
  normalizationScale = 1.0,
  targetUpdateFactor = 0,
  verbose = FALSE,
  maxit = 1e4
) {

  n <- nrow(xTrain)
  d <- ncol(xTrain)

  normalization <- PolyPropR::getNormalization(xTrain, normalizationType, normalizationScale)
  xTrain <- PolyPropR::normalize(xTrain, normalization)

  # Initial linear fit
  input1 <- xTrain[-n, , drop = FALSE]
  output1 <- xTrain[-1, , drop = FALSE]
  features1 <- PolyPropR::evaluateMonomialFeatures(input1, nDeg)
  coef <- PolyPropR::fitLinear(features1, output1)

  analysis <- xTrain
  target <- xTrain

  for (k in seq_along(weightSchedule)) {

    if (verbose) cat("Start weight", k, "\n")

    ws <- weightSchedule[[k]]

    obj <- TMB::MakeADFun(
      data = c(model = "poly_model_full",
               list(
                 obs = target,
                 weights_obs = ws$weightsObs,
                 weights_ana = ws$weightsAna,
                 weights_pen = ws$weightsPen,
                 weights_coef = ws$weightsCoef,
                 deg = nDeg,
                 intermediate = ws$intermediate
               )),
      parameters = list(coef = coef, ana = analysis),
      DLL = "DynLossOptimR_TMBExports",
      silent = TRUE
    )

    opt <- optimizer(obj$par, obj$fn, obj$gr)
    if (opt$convergence > 1L) break # 0: all good. 1: not converged but maybe good enough

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


#' @export
predictFullLoss <- function(model, initialConditions, nPred) {
  predictDefault(model, initialConditions, nPred)
}
