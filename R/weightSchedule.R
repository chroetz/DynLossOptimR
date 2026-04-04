getWeightScheduleFull <- function(
  nDeg,
  d,
  kmax = 1L,
  weightsObsBase = 0.2,
  weightsAnaBase = 0.2,
  weightsPen = double(0),
  reltol = 1e-7
) {
  weightSchedule <- lapply(seq_len(kmax), function(k) {
    list(
      weightsObs = c(1, rep(weightsObsBase, k)) / k,
      weightsAna = rep(weightsAnaBase * k, k) / k,
      weightsPen = weightsPen,
      weightsCoef = rep(0, nrow(PolyPropR::getMonomialFeatureDegrees(d, nDeg))),
      reltol = reltol
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


getWeightScheduleCoef <- function(
  nDeg,
  d,
  kmax = 1L,
  weightsObsBase = 0.2,
  reltol = 1e-7
) {
  weightSchedule <- lapply(seq_len(kmax), function(k) {
    list(
      weightsObs = c(1, rep(weightsObsBase, k)) / k,
      weightsCoef = rep(0, nrow(PolyPropR::getMonomialFeatureDegrees(d, nDeg))),
      reltol = reltol
    )
  })
  return(weightSchedule)
}

