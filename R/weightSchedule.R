#' @export
getWeightScheduleFull <- function(
  nDeg,
  d,
  intermediate = 1L,
  kmax = 1L,
  weightsObsBase = 0.2,
  weightsAnaBase = 0.2,
  weightsPen = double(0)
) {
  weightSchedule <- lapply(seq_len(kmax), function(k) {
    list(
      weightsObs = c(1, rep(weightsObsBase, k)) / k,
      weightsAna = rep(weightsAnaBase * k, k) / k,
      weightsPen = weightsPen,
      weightsCoef = rep(0, nrow(PolyPropR::getMonomialFeatureDegrees(d, nDeg))),
      intermediate = 1L
      # weightsCoef = dplyr::case_match(
      #   rowSums(PolyPropR::getMonomialFeatureDegrees(d, nDeg)),
      #   0 ~ 0,
      #   1 ~ 0,
      #   2 ~ 2e3,
      #   3 ~ 4e5
      # )
    )
  })
  return(weightSchedule)
}


#' @export
getWeightScheduleCoef <- function(
    nDeg,
    d,
    intermediate = 1L,
    kmax = 1L,
    weightsObsBase = 0.2,
    type = "const"
) {
  obsWeightFun <- switch(
    type,
    const = \(k, w) c(1, rep(w, k)) / k,
    point = \(k, w) c(rep(0, k), w),
    stop("Unknown weight schedule type ", type)
  )
  weightSchedule <- lapply(seq_len(kmax), function(k) {
    list(
      weightsObs = obsWeightFun(k, weightsObsBase),
      weightsCoef = rep(0, nrow(PolyPropR::getMonomialFeatureDegrees(d, nDeg))),
      intermediate = 1L
    )
  })
  return(weightSchedule)
}

#' @export
getWeightScheduleCoefIntermediate <- function(
    nDeg,
    d,
    intermediate,
    weightsObsBase = 1
) {
  weightSchedule <- lapply(seq_len(intermediate), function(k) {
    list(
      weightsObs = c(0, weightsObsBase),
      weightsCoef = rep(0, nrow(PolyPropR::getMonomialFeatureDegrees(d, nDeg))),
      intermediate = k
    )
  })
  return(weightSchedule)
}

#' @export
getWeightScheduleAna <- function(
    nDeg,
    d,
    kmax = 1L,
    weightsObsBase = 0.2,
    weightsAnaBase = 0.2,
    weightsPen = double(0)
) {
  weightSchedule <- lapply(seq_len(kmax), function(k) {
    list(
      weightsObs = c(1, rep(weightsObsBase, k)) / k,
      weightsAna = rep(weightsAnaBase * k, k) / k,
      weightsPen = weightsPen
    )
  })
  return(weightSchedule)
}

