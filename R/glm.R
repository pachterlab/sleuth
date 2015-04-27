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
          if (i %% 3000 == 0) {
            print(i)
          }
          structure(
            fit_gamma_glm(design_mat, exp_mat[i,]),
            class = c("glm", "lm")
            )
        },
        error = function(err) { err },
        finally = function() {}
        )
    }), rownames(exp_mat))
}

#' @export
bs_null_obs_glm_fit <- function(obj, downsample = NULL) {
  stopifnot( is(obj, "sleuth") )

  obs_norm <- spread_abundance_by(obj$obs_norm, "est_counts")
  bs <- dcast_bootstrap(obj, "est_counts", downsample)

  bs_design <- data.frame(bootstrap = c(rep.int( FALSE, ncol(obs_norm) ),
      rep.int(TRUE, ncol(bs))))
  all_data <- cbind(obs_norm, bs)
  rm(obs_norm, bs)

  design_mat <- model.matrix(~ bootstrap, bs_design)

  glm_by_rows(all_data, design_mat)
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

#' @export
fixed_null_glm_all_bs <- function(obj, downsample = 30) {
  stopifnot( is(obj, "sleuth") )

  obs_norm <- spread_abundance_by(obj$obs_norm, "est_counts")
  bs <- dcast_bootstrap(obj, "est_counts", downsample)

  model_mat <- model.matrix( ~ bio_rep,
    data = data.frame(bio_rep = gl(ncol(obs_norm), downsample)))


  glm_by_rows(bs, model_mat)
}


#' @export
fixed_null_glm_all_bs_no_grouping <- function(obj, downsample = 30) {
  stopifnot( is(obj, "sleuth") )

  obs_norm <- spread_abundance_by(obj$obs_norm, "est_counts")
  bs <- dcast_bootstrap(obj, "est_counts", downsample)

  model_mat <- model.matrix(~1,
    data = data.frame(bio_rep = gl(ncol(obs_norm), downsample)))

  glm_by_rows(bs, model_mat)
}
