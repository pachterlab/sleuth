#' @export
fit_gamma_glm <- function(design_mat, response) {
  # TODO: incorporate offset
  # TODO: figure out how to get dispersion param back
  glm.fit(x = design_mat, y = response, family = Gamma(link = "log"))
}


#' @export
glm_by_rows <- function(exp_mat, design_mat) {
  setNames(lapply(1:nrow(exp_mat),
    function (i)
    {
      tryCatch (
        {
          if (i %% 1000 == 0) {
            print(i)
          }
          structure(
            fit_gamma_glm(design_mat, exp_mat[i,]),
            class = c("glm", "lm")
            )
        },
        error = function(err) {},
        finally = function() {}
        )
    }), rownames(exp_mat))
}

#' @export
bs_obs_glm_fit <- function(obj) {
  stopifnot( is(obj, "sleuth") )
}

#' @export
fit_null_glms <- function(obj) {
  stopifnot( is(obj, "sleuth") )

  obs_norm <- spread_abundance_by(obj$obs_norm, "est_counts")

  cat("Fitting null GLMs\n")
  # fit the intercept only model
  obj$valid_null_glms <- glm_by_rows(obs_norm,
    model.matrix(~ 1, s_o$sample_to_condition)) %>%
      Filter(function(x) is(x, "glm"), .)

  cat("Computing summaries")
  obj$null_glm_summary <- Map(summary, obj$valid_null_glms)

  obj
}

#' @export
locfit_phi <- function(obj) {
  stopifnot( is(obj, "sleuth") )
  stopifnot( !is.null(obj[["null_glm_summary"]]) )

  disp <- unlist(Map(function(x) x[["dispersion"]], obj$null_glm_summary))
  obs_norm <- spread_abundance_by(obj$obs_norm, "est_counts")
  null_mean <- apply(obs_norm, 1, mean)
  null_mean <- null_mean[names(disp)]

  list(null_mean, disp)
}

#' @export
glm_coef_by_rows <- function(exp_mat, design_mat) {
  res <- matrix(NA_real_,
    nrow = nrow(exp_mat), ncol = ncol(design_mat),
    dimnames = list(rownames(exp_mat), colnames(design_mat)))

  for (i in 1:nrow(exp_mat)) {
    tryCatch (
      {
        res[i,] <- coef(fit_gamma_glm(design_mat, exp_mat[i, ]))
      },
      error = function(err) {},
      finally = function() {}
      )
  }

  res
}

merge_glms <- function(row_glm) {
  stopifnot( is(row_glm, "list") )

  n_trans <- length(row_glm[[1]])


}
