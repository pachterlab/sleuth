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
print.kallisto <- function(obj) {
  cat("\tkallisto object\n")
  cat("\n")
  cat("transcripts: ", length(obj$abundance$target_id), "\n")
  n_bootstrap <- length(obj$bootstrap)

  cat("bootstraps: ", n_bootstrap, "\n")

  invisible(obj)
}

#' @export
head.kallisto <- function(obj) {
  print(obj)
}
