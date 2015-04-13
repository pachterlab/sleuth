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
