optimizer <- function(par, fn, gr, maxit = 100, factr = 1e7, pgtol = 1e-5) {
  result <- tryCatch(
    stats::optim(par, fn, gr,
                 method = "L-BFGS-B",
                 control = list(
                   maxit = maxit,
                   factr = factr,
                   pgtol = pgtol
                 )
    ),
    error = function(e) {
      list(
        par         = par,
        value       = NA_real_,
        counts      = c("function" = NA_integer_, gradient = NA_integer_),
        convergence = 99L,
        message     = conditionMessage(e)
      )
    }
  )

  if (result$convergence != 0L) {
    warning("Optimizer did not converge (code ", result$convergence, "): ",
            result$message, immediate. = TRUE)
  }

  result
}
