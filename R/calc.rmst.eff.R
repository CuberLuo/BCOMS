source("get.post.sample.R")
source("lomax.distr.R")

calc.rmst.eff <- function(pfs.data, immune.data, dose.data, JJ, KK, tau, silent = TRUE) {
  ndose <- JJ * KK
  dose.jags <- get.post.sample(
    pfs.data = pfs.data,
    immune.data = immune.data,
    dose.data = dose.data,
    n.pat = length(pfs.data),
    ndose = ndose,
    tau = tau,
    silent = silent
  )
  rmst.matrix <- matrix(0, JJ, KK)
  pE.matrix <- matrix(0, JJ, KK)

  for (d in seq_len(ndose)) {
    alpha.post <- mean(dose.jags[, paste0("alpha[", d, "]")])
    gamma.post <- mean(dose.jags[, paste0("gamma[", d, "]")])
    pE.post <- mean(dose.jags[, paste0("pE[", d, "]")])

    rmst.matrix[d] <- integrate(slomax,
      lower = 0, upper = tau,
      alpha = alpha.post, gamma = gamma.post
    )$value
    pE.matrix[d] <- pE.post
  }
  return(list(rmst = rmst.matrix, pE = pE.matrix))
}
