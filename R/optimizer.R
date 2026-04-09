optimizer <- function(par, fn, gr, he = NULL, method = c("nlminb", "L-BFGS-B"),
                      maxit = 1e4, restarts = 2) {
  method <- match.arg(method)

  run <- function(par) {
    if (method == "nlminb") {
      res <- nlminb(par, fn, gr, hessian = he,
                    control = list(iter.max = maxit, eval.max = 2 * maxit,
                                   rel.tol = 1e-10))
      list(
        par         = res$par,
        value       = res$objective,
        counts      = c("function" = res$evaluations[["function"]],
                        gradient  = res$evaluations[["gradient"]]),
        convergence = res$convergence,
        message     = res$message
      )
    } else {
      stats::optim(par, fn, gr,
                   method = "L-BFGS-B",
                   control = list(maxit = maxit, factr = 1e6, pgtol = 1e-6))
    }
  }

  result <- NULL
  for (i in seq_len(restarts + 1)) {
    result <- tryCatch(run(par), error = function(e) {
      list(
        par         = par,
        value       = NA_real_,
        counts      = c("function" = NA_integer_, gradient = NA_integer_),
        convergence = 99L,
        message     = conditionMessage(e)
      )
    })
    if (result$convergence == 0L) break
    par <- result$par
  }

  if (result$convergence != 0L) {
    warning("Optimizer did not converge (code ", result$convergence, "): ",
            result$message, immediate. = TRUE)
  }
  result
}
