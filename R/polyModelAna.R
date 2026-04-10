#' @export
fitAnaLoss <- function(
  xTrain,
  nDeg,
  coef,
  weightSchedule,
  intermediate = 1L,
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

  analysis <- xTrain
  target <- xTrain

  for (k in seq_along(weightSchedule)) {

    if (verbose) cat("Start weight", k, "\n")

    ws <- weightSchedule[[k]]

    obj <- TMB::MakeADFun(
      data = c(model = "poly_model_ana",
               list(
                 obs = target,
                 coef = coef,
                 weights_obs = ws$weightsObs,
                 weights_ana = ws$weightsAna,
                 weights_pen = ws$weightsPen,
                 deg = nDeg,
                 intermediate = intermediate
               )),
      parameters = list(ana = analysis),
      DLL = "DynLossOptimR_TMBExports",
      silent = TRUE
    )

    opt <- optimizer(obj$par, obj$fn, obj$gr)
    if (opt$convergence > 1L) break # 0: all good. 1: not converged but maybe good enough

    analysisNew <- opt$par[names(opt$par) == "ana"]
    dim(analysisNew) <- dim(analysis)

    analysis <- analysisNew
    target <- (1 - targetUpdateFactor) * target + targetUpdateFactor * analysis
  }

  return(
    list(
      coef = coef,
      nDeg = nDeg,
      analysis = PolyPropR::denormalize(analysis, normalization),
      normalization = normalization,
      intermediate = intermediate
    )
  )

}


#' @export
predictAnaLoss <- function(model, initialConditions, nPred) {
  predictDefault(model, initialConditions, nPred)
}
