source("select.mtdc.subtrial.R")
source("lomax.distr.R")
comb.subtrial <- function(target, pT.true, pE.true, alpha.true, gamma.true, dosespace,
                          npts, ntox, neff, elimi, ncohort, cohortsize, subtriali,
                          n.earlystop, startdose, cutoff.eli, p.saf, p.tox,
                          maxtT = 1, maxtE = 1, accrual = 3, arrival.dist = 1) {
  pT.truee <- pT.true[dosespace]
  pE.truee <- pE.true[dosespace]
  gamma.truee <- gamma.true[dosespace]
  ndose <- length(dosespace)
  selectdose <- 0
  lambda1 <- log((1 - p.saf) / (1 - target)) / log(target * (1 - p.saf) / (p.saf * (1 - target)))
  lambda2 <- log((1 - target) / (1 - p.tox)) / log(p.tox * (1 - target) / (target * (1 - p.tox)))
  yT <- yE <- n <- rep(0, ndose)
  pfs.data <- NULL
  immune.data <- NULL
  dose.data <- NULL
  earlystop <- 0
  d <- startdose
  elmT <- rep(0, ndose)
  t.enter <- NULL
  t.decision <- 0
  duration <- 0
  # cat("-------------------------\n")
  # cat("subtriali:", subtriali, "\n")
  for (icohort in 1:ncohort) {
    for (j in 1:cohortsize) {
      if (j == 1) {
        t.enter <- c(t.enter, t.decision)
      } else {
        if (arrival.dist == 1) {
          t.enter <- c(t.enter, t.enter[length(t.enter)] + runif(1, 0, 2 / accrual))
        } else if (arrival.dist == 2) {
          t.enter <- c(t.enter, t.enter[length(t.enter)] + rexp(1, accrual))
        }
      }
    }
    yT.outcome <- rbinom(cohortsize, 1, pT.truee[d])
    yT[d] <- yT[d] + sum(yT.outcome)
    # cat("dose", d, "yT.outcome", yT.outcome, "\n")
    yE.outcome <- rbinom(cohortsize, 1, pE.truee[d])
    yE[d] <- yE[d] + sum(yE.outcome)
    pfs.data <- c(pfs.data, rlomax(cohortsize, alpha.true, gamma.truee[d]))
    immune.data <- c(immune.data, yE.outcome)
    dose.data <- c(dose.data, rep(d, cohortsize))
    n[d] <- n[d] + cohortsize
    pT.dose <- yT[d] / n[d]

    if (icohort == ncohort) {
      t.decision <- t.enter[length(t.enter)]
    } else {
      t.decision <- t.enter[length(t.enter)] + max(maxtT, maxtE)
    }

    # 当前剂量毒性过高，淘汰当前剂量以及更高剂量
    if (1 - pbeta(target, yT[d] + 1, n[d] - yT[d] + 1) > cutoff.eli) {
      elmT[d:ndose] <- 1
      if (d == 1) {
        earlystop <- 1
        break
      }
    }

    if (pT.dose <= lambda1 && d != ndose) {
      if (elmT[d + 1] == 0) {
        d <- d + 1
      }
    } else if (pT.dose >= lambda2 && d != 1) {
      d <- d - 1
    } else {
      d <- d
    }
    if (n[d] >= n.earlystop) break
  }
  duration <- t.decision

  if (earlystop == 1) {
    selectdose <- 99
  } else {
    # 根据n,yT选择当前子试验的候选MTDC
    selectdose <- select.mtdc.subtrial(
      target = target, n = n, yT = yT, elmT = elmT
    )
  }
  # cat("selectdose:", selectdose, "\n")
  npts[dosespace] <- n
  ntox[dosespace] <- yT
  neff[dosespace] <- yE
  elimi[dosespace] <- elmT
  list(
    npts = npts, ntox = ntox, neff = neff, elimi = elimi,
    totaln = sum(npts), selectdose = selectdose, duration = duration,
    pfs.data = pfs.data, immune.data = immune.data, dose.data = dosespace[dose.data]
  )
}
