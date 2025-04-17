generate_split_dataset <- function(data, size) {

  if (size < 0 || size > 1) {
    stop("Size must be a number between 0 and 1.")
  }

  data |>
    dplyr::sample_frac(size) |>
    dplyr::pull(rowId) |>
    sort()

}


