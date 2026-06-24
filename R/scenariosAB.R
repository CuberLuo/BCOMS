pT.trueA <- matrix(c(
  0.35, 0.38, 0.40,
  0.37, 0.41, 0.42
), nrow = 2, byrow = TRUE)

pE.trueA <- matrix(c(
  0.10, 0.15, 0.18,
  0.14, 0.17, 0.25
), nrow = 2, byrow = TRUE)

rmst.trueA <- matrix(c(
  2.12, 2.31, 2.40,
  2.29, 2.27, 2.76
), nrow = 2, byrow = TRUE)

pT.trueB <- matrix(c(
  0.01, 0.03, 0.05,
  0.03, 0.06, 0.07
), nrow = 2, byrow = TRUE)

pE.trueB <- matrix(c(
  0.15, 0.24, 0.26,
  0.26, 0.28, 0.52
), nrow = 2, byrow = TRUE)

rmst.trueB <- matrix(c(
  4.40, 5.04, 5.08,
  5.20, 5.18, 7.07
), nrow = 2, byrow = TRUE)


list.p.true.tox <- list(
  pT.trueA, pT.trueB
)
list.p.true.eff <- list(
  pE.trueA, pE.trueB
)
list.rmst.true <- list(
  rmst.trueA, rmst.trueB
)
list.tdose <- list(
  c(-1), c(6)
)
cat("Loaded scenarios!\n")
