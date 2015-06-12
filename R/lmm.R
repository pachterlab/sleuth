#' @export
fit_lmm <- function(design_mat, design_formula, random_formula, result) {
  stopifnot( length(result) == nrow(design_mat) )

  design_mat <- design_mat %>%
    mutate(expression = result)

  nlme::lme(design_formula, random = random_formula, data = design_mat)
}

#' @export
lmm_design <- function(obj) {
  stopifnot( is(obj, "sleuth") )

  # rownames on DFs are annoying...
  result <- tbl_df(obj$sample_to_condition)

  n_bs <- sapply(obj$kal, function(x) length(x$bootstrap))

  bs_ids <- rep.int(1:nrow(result), n_bs)

  result[bs_ids,]
}

#' @export
lmm_by_row <- function(obj, design_formula, random_formula, filter_df) {
  stopifnot( is(obj, "sleuth") )

  filter_df <- as.data.frame(filter_df)
  rownames(filter_df) <- filter_df$target_id
  filter_df$target_id <- NULL

  bs <- dcast_bootstrap(obj, "est_counts")

  design <- lmm_design(obj)

  targs <- rownames(bs)
  filt <- filter_df[targs,1]
  stopifnot(length(targs) == length(filt))

  bs_lmm <- lapply(1:nrow(bs),
    function(i)
    {
      tryCatch(
        {
          if (filt[i]) {
            fit_lmm(design, design_formula, random_formula, bs[i,])
          } else {
            stop("didn't pass filter")
          }
        },
        error = function(e) {e},
        finally = function() {}
        )
    }
    )

  bs_lmm
}

#' @export
summary_lme <- function (object, stdFixed, verbose = FALSE, ...) {
  fixed <- fixef(object)
  # stdFixed <- sqrt(diag(as.matrix(object$varFix)))
  # object$corFixed <- array(t(object$varFix/stdFixed)/stdFixed,
  #     dim(object$varFix), list(names(fixed), names(fixed)))
  # if (adjustSigma && object$method == "ML")
  #     stdFixed <- stdFixed * sqrt(object$dims$N/(object$dims$N -
  #         length(stdFixed)))

  tTable <- data.frame(fixed, stdFixed, object$fixDF[["X"]],
    fixed/stdFixed, fixed)
  dimnames(tTable) <- list(names(fixed), c("Value", "Std.Error",
      "DF", "t-value", "p-value"))
  tTable[, "p-value"] <- 2 * pt(-abs(tTable[, "t-value"]),
    tTable[, "DF"])
  object$tTable <- as.matrix(tTable)
  resd <- resid(object, type = "pearson")
  if (length(resd) > 5) {
    resd <- quantile(resd, na.rm = TRUE)
    names(resd) <- c("Min", "Q1", "Med", "Q3", "Max")
  }
  object$residuals <- resd
  aux <- logLik(object)
  object$BIC <- BIC(aux)
  object$AIC <- AIC(aux)
  attr(object, "oClass") <- class(object)
  attr(object, "verbose") <- verbose
  class(object) <- c("summary.lme", class(object))
  object
}
