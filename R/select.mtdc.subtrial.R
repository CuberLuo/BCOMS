select.mtdc.subtrial <- function(target, n, yT, elmT) {
  # cat("n:", n, "\n")
  # cat("yT:", yT, "\n")
  # cat("elmT:", elmT, "\n")
  elimi <- elmT
  if (elmT[1] == 1 || sum(n[elmT == 0]) == 0) {
    selectdose <- 99
  } else {
    adm.set <- (n != 0) & (elimi == 0)
    adm.index <- which(adm.set == TRUE)
    yT.adm <- yT[adm.set]
    n.adm <- n[adm.set]

    if (length(n.adm) == 0) {
      selectdose <- 99
    } else {
      pT <- (yT.adm + 0.05) / (n.adm + 0.1) # II期毒性先验Beta(0.05,0.05)
      pT <- pava(pT, n.adm + 0.1) + 0.001 * seq_len(length(pT))
      # cat("pT:", round(pT, 2), "\n")
      d_mtd.selected <- which.min(abs(pT - target)) # MTD
      selectdose <- adm.index[d_mtd.selected]
    }
  }

  return(selectdose)
}
