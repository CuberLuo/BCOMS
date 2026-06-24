source("select.mtdc.subtrial.R")
source("gen.tite.R")
source("lomax.distr.R")

comb.subtrial.tite <- function(target, pT.true, pE.true, alpha.true, gamma.true, dosespace,
                               npts, ntox, neff, elimi, ncohort, cohortsize, subtriali,
                               n.earlystop, startdose, cutoff.eli, p.saf, p.tox, distT = 2,
                               distE = 2, maxtT = 1, maxtE = 1, accrual = 3, arrival.dist = 1) {
  pT.truee <- pT.true[dosespace]
  pE.truee <- pE.true[dosespace]
  gamma.truee <- gamma.true[dosespace]
  ndose <- length(dosespace)
  selectdose <- 0
  lambda1 <- log((1 - p.saf) / (1 - target)) / log(target * (1 - p.saf) / (p.saf * (1 - target)))
  lambda2 <- log((1 - target) / (1 - p.tox)) / log(p.tox * (1 - target) / (target * (1 - p.tox)))
  yT <- yE <- n <- rep(0, ndose)
  nT.ess <- rep(0, ndose)
  pfs.data <- NULL
  immune.data <- NULL
  dose.data <- NULL
  earlystop <- 0
  d <- startdose
  elmT <- rep(0, ndose)
  tT.enter <- NULL
  tT.event <- NULL
  tT.decision <- 0
  yT.i <- yE.i <- NULL
  duration <- 0
  # cat("\n subtriali:", subtriali, "\n")

  for (icohort in 1:ncohort) {
    # cat("----------------------------------\n")
    # cat("Cohort:", icohort, "Dose:", d, "\n")

    for (j in 1:cohortsize) {
      if (j == 1) {
        tT.enter <- c(tT.enter, tT.decision)
      } else {
        if (arrival.dist == 1) {
          tT.enter <- c(tT.enter, tT.enter[length(tT.enter)] + runif(1, 0, 2 / accrual))
        } else if (arrival.dist == 2) {
          tT.enter <- c(tT.enter, tT.enter[length(tT.enter)] + rexp(1, accrual))
        }
      }
    }

    # cat("tT.enter:", tail(tT.enter*30, cohortsize), "\n")

    obs.tox <- gen.tite(cohortsize, pT.truee[d], distT, Tobs = maxtT)
    obs.eff <- gen.tite(cohortsize, pE.truee[d], distE, Tobs = maxtE)
    yT.outcome <- obs.tox$event
    yE.outcome <- obs.eff$event
    # cat("yT.outcome:", yT.outcome, "yE.outcome:", yE.outcome, "\n")

    pfs.data <- c(pfs.data, rlomax(cohortsize, alpha.true, gamma.truee[d]))
    immune.data <- c(immune.data, yE.outcome)
    dose.data <- c(dose.data, rep(d, cohortsize))
    # 实际的毒性和疗效结果
    yT.i <- c(yT.i, yT.outcome)
    yE.i <- c(yE.i, yE.outcome)
    # 毒性和疗效的event time
    tT.event <- c(tT.event, obs.tox$t.event)
    # tE.event <- c(tE.event, obs.eff$t.event)
    # cat("tT.event:", tail(tT.event, cohortsize), "\n")
    # 初始化当前cohort的decision time
    tT.decision <- tT.enter[length(tT.enter)]
    pending <- 1

    while (pending == 1) {
      pending <- 0
      # 下一个cohort患者到来的时刻即为当前cohort的decision time
      if (icohort == ncohort) {
        tT.decision <- tT.decision + maxtT
      } else {
        if (arrival.dist == 1) {
          tT.decision <- tT.decision + runif(1, 0, 2 / accrual)
        } else if (arrival.dist == 2) {
          tT.decision <- tT.decision + rexp(1, accrual)
        }
      }

      # cat("decision time:", tT.decision*30, "\n")

      deltaT <- (tT.enter + tT.event <= tT.decision) # true表示观察到毒性了
      # cat("deltaT:", tail(deltaT, cohortsize), "\n")
      # cat("update nT.ess for each dose:\n")

      for (dd in 1:ndose) {
        idx <- which(dose.data == dd)
        if (length(idx) == 0) next
        yT[dd] <- sum(yT.i[idx][deltaT[idx]])
        yE[dd] <- sum(yE.i[idx])
        pendingT.weight <- sum(
          pmin(tT.decision - tT.enter[idx][!deltaT[idx]], maxtT)
        ) / maxtT
        nT.ess[dd] <- sum(deltaT[idx]) + pendingT.weight
        n[dd] <- length(idx)
      }

      # 至少一半患者完成毒性
      if (sum(deltaT[dose.data == d]) >= 0.5 * n[d]) {
        # nc <- ceiling(n[d] / cohortsize)
        # nc <- max(1, ceiling(nT.ess[d] / cohortsize))
        pT.dose <- yT[d] / nT.ess[d]
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
      } else {
        pending <- 1
      }
    }

    if (n[d] >= n.earlystop) break
  }

  duration <- tT.decision
  # cat("duration:", duration*30, "\n")

  if (earlystop == 1) {
    selectdose <- 99
  } else {
    # 根据n,yT选择当前子试验的候选MTDC
    selectdose <- select.mtdc.subtrial(
      target = target, n = n, yT = yT, elmT = elmT
    )
  }

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
