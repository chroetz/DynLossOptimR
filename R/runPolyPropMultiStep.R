#' Multi-Step Polynomial Propagation
#'
#' Fits a polynomial map by progressively optimizing over longer forward
#' prediction horizons using TMB for automatic differentiation, then
#' generates a forecast by iterating the learned map.
#'
#' @param xTrain Numeric matrix (n x d) of training data, where each row
#'   is a state vector.
#' @param nPred Integer; number of forecast steps to produce.
#' @param nDeg Integer; degree of polynomial features (default 3).
#' @param forwardSteps Integer vector of forward-step horizons to optimize
#'   over sequentially (default \code{c(2, 4, 8)}).
#'
#' @return A list with components:
#'   \describe{
#'     \item{assimilation}{The training data (unchanged).}
#'     \item{forecast}{Numeric matrix (nPred x d) of predicted states.}
#'   }
#'
#' @export
runPolyPropMultiStep <- function(xTrain, nPred,
                                  nDeg = 3L,
                                  forwardSteps = c(2, 4, 8)) {

  n <- nrow(xTrain)
  d <- ncol(xTrain)

  optimizer <- function(par, fn, gr) {
    stats::optim(par, fn, gr,
      method = "BFGS",
      control = list(
        maxit = 1e5,
        reltol = 1e-14
      )
    )
  }

  # Initial linear fit
  input1 <- xTrain[-n, , drop = FALSE]
  output1 <- xTrain[-1, , drop = FALSE]
  features1 <- PolyPropR::evaluateMonomialFeatures(input1, nDeg)
  coef <- PolyPropR::fitLinear(features1, output1)

  # Progressive multi-step refinement via TMB
  for (k in forwardSteps) {
    input <- xTrain[seq_len(n - k), , drop = FALSE]
    output <- xTrain[(k + 1):n, , drop = FALSE]

    obj <- TMB::MakeADFun(
      data = c(model = "poly_model",
               list(X = input, Y = output, deg = nDeg, weights = rep(1, k))),
      parameters = list(Theta = coef),
      DLL = "DynLossOptimR_TMBExports",
      silent = TRUE
    )

    opt <- tryCatch(
      optimizer(obj$par, obj$fn, obj$gr),
      error = function(cond) NULL
    )

    if (is.null(opt)) {
      message("Optimization failed at forward step ", k)
      break
    }
    coef <- matrix(opt$par, ncol = d)
  }

  # Generate forecast
  prediction <- matrix(NA_real_, nrow = nPred, ncol = d)
  currentState <- xTrain[n, ]
  for (i in seq_len(nPred)) {
    currentState <- PolyPropR::evaluateMonomialFeatures(matrix(currentState, nrow = 1), nDeg) %*% coef
    prediction[i, ] <- currentState
  }

  list(
    assimilation = xTrain,
    forecast = prediction
  )
}
