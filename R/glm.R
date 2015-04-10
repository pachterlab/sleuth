#' @export
fit_gamma_glm <- function(design_mat, response) {
  # TODO: incorporate offset
  # TODO: figure out how to get dispersion param back
  glm.fit(x = design_mat, y = response, family = Gamma(link = "log"))
}
