#' Get labfs from betas and standard errors
#'
#' @param effect.est matrix of snp regression coefficients (i.e. regression beta
#'  values) in the genomic region
#' @param effect.se matrix of standard errors associated with the beta values
#' @param binary.outcomes a binary vector of dimension the number of traits:
#'  1 represents a binary trait 0 otherwise
#' @return A matrix containing labfs
#' @export
#'
pg_get_labfs <- function(effect.est, effect.se, binary.outcomes) {
  if (any(is.na(effect.est))) {
    stop("there are missing values in effect.est")
  }
  if (any(is.na(effect.se))) {
    stop("there are missing values in effect.se")
  }
  if (0 %in% effect.se) {
    stop("there are zero values in effect.se")
  }
  if(!is.vector(binary.outcomes)) {
    # NOTE: In R, both 1 and c(1) are vectors!  Therefore, by doing is.vector,
    # we are making sure that binary.outcomes is either a single value or a
    # vector of values.
    stop("binary.outcomes should be a vector or a single value")
  }
  if (any(is.na(binary.outcomes))) {
    stop("there are missing values in binary.outcomes")
  }
  if (any(!(binary.outcomes %in% c(0, 1)))) {
    stop("binay.outcomes should contain values that are either 0s or 1s")
  }
  if(!(is.matrix(effect.est) == is.matrix(effect.se))) {
    stop("effect.est and effect.se should both be either matrices or vectors,
          not a mix of the two.")
  }
  if(!(is.atomic(effect.est) & is.atomic(effect.se))) {
    stop("effect.est and effect.se should both be either matrices or vectors")
  } else {
    if(!(is.matrix(effect.est) == is.matrix(effect.se))) {
      stop("effect.est and effect.se should both be either matrices or vectors,
          not a mix of the two.")
    }
  }
  if(is.matrix(effect.est)) {
    if(length(binary.outcomes) != dim(effect.est)[[2]]) {
      stop("as effect.est and effect.st are matrices, binary.outcomes should be
        a vector with it's length matching the number of columns in the
        matrices")
    }
  } else {
    if(length(binary.outcomes) != 1) {
      stop("as effect.est and effect.st are vectors, binary.outcomes should be
        a single value or a vector of length 1")
    }
  }

  Z <- effect.est / effect.se

  W <- 0.15 / effect.se

  if(is.matrix(effect.se)) {
    W[, which(binary.outcomes == 1)] <- 0.2 / (effect.se[, which(binary.outcomes == 1)])
  } else {
    if(binary.outcomes == 1) {
      W <- 0.2 / (effect.se)
    }
  }

  Zsq <- Z^2
  Wsq <- 1 / (1 + W^2)

  labf.result <- pg_rapid_labf(Zsq, Wsq)
  return(labf.result)
}

#' Rapidly calculate log approximate Bayes Factors (labf)
#'
#' @param Zsq matrix of Z-scores
#' @param Wsq ratio matrix of prior standard deviation and observed standard
#'  errors squared
#' @return matrix of labfs
#' @export
pg_rapid_labf <- function(Zsq, Wsq) {
    return(
        0.5 * (log(Wsq) + (Zsq) * (1 - Wsq))
    )
}