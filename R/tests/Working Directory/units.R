#
#    sleuth: inspect your RNA-Seq with a pack of kallistos
#
#    Copyright (C) 2015  Harold Pimentel, Nicolas Bray, Pall Melsted, Lior Pachter
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' @export
counts_to_tpm <- function(est_counts, eff_len) {
  stopifnot( length(eff_len) == length(est_counts) )

  which_valid <- which(eff_len > 0)

  num <- (est_counts / eff_len)
  num[-which_valid] <- 0
  denom <- sum(num)

  (1e6 * num) / denom
}

#' @export
tpm_to_alpha <- function(tpm, eff_len) {
  stopifnot( length(tpm) == length(eff_len) )

  num <- tpm * eff_len
  denom <- sum(num)

  num / denom
}

#' @export
counts_to_fpkm <- function(counts, eff_len) {
  stopifnot( length(counts) == length(eff_len) )

  which_valid <- which(eff_len > 0)
  N <- sum(counts)
  fpkm <- (counts / eff_len) * (1e9 / N)
  fpkm[-which_valid] <- 0

  fpkm
}
