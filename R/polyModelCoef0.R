#' @export
fitCoef0Loss <- function(
  xTrain,
  nDeg,
  layers,
  normalizationType = "none",
  normalizationScale = 1.0,
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

  for (k in seq_len(layers-1)) {

    if (verbose) cat("Start ", k, "\n")

    obj <- TMB::MakeADFun(
      data = c(model = "poly_model_coef0",
               list(
                 obs = target,
                 deg = nDeg,
                 layers = k+1
               )),
      parameters = list(coef = coef),
      DLL = "DynLossOptimR_TMBExports",
      silent = TRUE
    )

    opt <- optimizer(obj$par, obj$fn, obj$gr)
    if (opt$convergence > 1L) break # 0: all good. 1: not converged but maybe good enough

    coefNew <- opt$par[names(opt$par) == "coef"]
    dim(coefNew) <- dim(coef)

    coef <- coefNew
  }

  return(
    list(
      coef = coef,
      nDeg = nDeg,
      analysis = PolyPropR::denormalize(xTrain, normalization),
      normalization = normalization,
      intermediate = 1
    )
  )

}


#' @export
predictCoef0Loss <- function(model, initialConditions, nPred) {
  predictDefault(model, initialConditions, nPred)
}
