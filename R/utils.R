one2two <- function(index, row.num, col.num) {
  col.index <- ifelse(index %% row.num == 0, index %/% row.num,
    index %/% row.num + 1
  )
  row.index <- index - (col.index - 1) * row.num
  return(c(row.index, col.index))
}

two2one <- function(row.index, col.index, row.num) {
  return((col.index - 1) * row.num + row.index)
}

format.matrix <- function(mat, digits = 2) {
  formatted_str <- formatC(mat, digits = digits, format = "f", width = 5)
  numeric_vec <- as.numeric(formatted_str)
  result_mat <- matrix(numeric_vec, nrow = nrow(mat), ncol = ncol(mat))
  return(result_mat)
}

format.num <- function(number, digits = 1) {
  formatted_str <- formatC(number, digits = 1, format = "f")
  numeric_vec <- as.numeric(formatted_str)
  return(numeric_vec)
}

format.percent <- function(value, total) {
  percent <- 100 * value / total
  sprintf("%.1f%%", percent)
}

convert_to_matrix <- function(x) {
  if (is.numeric(x)) {
    if (length(x) == 1) {
      x <- matrix(x, nrow = 1, ncol = 1)
    } else if (!is.matrix(x)) {
      x <- matrix(x, ncol = 1)
    }
  }
  return(x)
}
