#' Full-Loss Polynomial Optimization with Analysis Update
#'
#' Jointly optimizes polynomial map coefficients and an analysis trajectory
#' by minimizing a composite loss that balances observation fit, analysis
#' smoothness, and coefficient regularization. Uses TMB for automatic
#' differentiation.
#'
#' @param xTrain Numeric matrix (n x d) of training observations.
#' @param nPred Integer; number of forecast steps to produce.
#' @param nDeg Integer; degree of polynomial features (default 3).
#' @param kmax Integer; number of optimization iterations (default 10).
#' @param targetUpdateFactor Numeric in \[0,1\]; blending factor for
#'   updating the target trajectory between iterations (default 0).
#' @param weightsObsBase Numeric; base weight for observation loss terms
#'   (default 0.2).
#' @param weightsAnaBase Numeric; base weight for analysis loss terms
#'   (default 0.2).
#' @param weightsPen Numeric vector; penalty weights (default
#'   \code{double(0)}).
#' @param reltol Numeric; relative tolerance for the optimizer
#'   (default 1e-12).
#'
#' @return A list with components:
#'   \describe{
#'     \item{assimilation}{Numeric matrix (n x d) of optimized analysis
#'       states.}
#'     \item{forecast}{Numeric matrix (nPred x d) of predicted states.}
#'   }
#'
#' @export
runFullLoss <- function(xTrain, nPred,
                        nDeg = 3L,
                        kmax = 10L,
                        targetUpdateFactor = 0,
                        weightsObsBase = 0.2,
                        weightsAnaBase = 0.2,
                        weightsPen = double(0),
                        reltol = 1e-12) {

  n <- nrow(xTrain)
  d <- ncol(xTrain)

  # Build per-iteration weight schedules
  weightSchedule <- lapply(seq_len(kmax), function(k) {
    list(
      weightsObs = c(1, rep(weightsObsBase, k)) / k,
      weightsAna = rep(weightsAnaBase * k, k) / k,
      reltol = reltol,
      weightsCoef = dplyr::case_match(
        rowSums(PolyPropR::getMonomialFeatureDegrees(d, nDeg)),
        0 ~ 0,
        1 ~ 0,
        2 ~ 2e3,
        3 ~ 4e5
      )
    )
  })

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

  for (k in seq_len(kmax)) {
    ws <- weightSchedule[[k]]

    obj <- TMB::MakeADFun(
      data = c(model = "poly_model_full",
               list(
                 obs = target,
                 weights_obs = ws$weightsObs,
                 weights_ana = ws$weightsAna,
                 weights_pen = weightsPen,
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

  # Generate forecast
  prediction <- matrix(NA_real_, nrow = nPred, ncol = d)
  currentState <- analysis[n, ]
  for (i in seq_len(nPred)) {
    currentState <- PolyPropR::evaluateMonomialFeatures(matrix(currentState, nrow = 1), nDeg) %*% coef
    prediction[i, ] <- currentState
  }

  list(
    assimilation = analysis,
    forecast = prediction
  )
}
