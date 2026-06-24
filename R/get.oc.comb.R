library(Iso)
library(parallel)
library(future.apply)
library(progressr)
library(nleqslv)
library(rjags)
library(runjags)
library(coda)
source("comb.subtrial.R")
source("comb.subtrial.tite.R")
source("utils.R")
source("calc.rmst.eff.R")

get.oc.comb <- function(pT.true, pE.true, rmst.true, alpha.true, gamma.true, target,
                        ncohort, cohortsize, n.earlystop, startdose = c(1, 1),
                        ntrial = 1000, cutoff.eli = 0.95, p.saf = 0.75 * target,
                        p.tox = 1.25 * target, lower_e, rmst.min, dose.max.sample,
                        tau, tite = FALSE) {
  JJ <- nrow(pT.true)
  KK <- ncol(pT.true)
  ntrial.obdc <- vector("list", ntrial)
  ntrial.nt <- vector("list", ntrial)
  ntrial.yt <- vector("list", ntrial)
  ntrial.ye <- vector("list", ntrial)
  ntrial.duration <- numeric(ntrial)
  init.startdose <- startdose

  run_subtrial <- function(dosespace, startdose, subtriali, state) {
    subtrial.fn <- if (tite) comb.subtrial.tite else comb.subtrial
    subtrial <- subtrial.fn(
      target = target, pT.true = pT.true, pE.true = pE.true,
      alpha.true = alpha.true, gamma.true = gamma.true, dosespace = dosespace,
      npts = state$npts, ntox = state$ntox, neff = state$neff, elimi = state$elimi,
      ncohort = ncohort[subtriali], cohortsize = cohortsize, subtriali = subtriali,
      n.earlystop = n.earlystop, startdose = startdose,
      cutoff.eli = cutoff.eli, p.saf = p.saf, p.tox = p.tox
    )

    list(
      subtrial = subtrial,
      selectdose = ifelse(subtrial$selectdose == 99, 99,
        dosespace[subtrial$selectdose]
      )
    )
  }

  update_state <- function(state, subtrial) {
    state$totaln <- subtrial$totaln
    state$npts <- subtrial$npts
    state$ntox <- subtrial$ntox
    state$neff <- subtrial$neff
    state$elimi <- subtrial$elimi
    state$total.duration <- state$total.duration + subtrial$duration
    state$stageI.pfs <- c(state$stageI.pfs, subtrial$pfs.data)
    state$stageI.immune <- c(state$stageI.immune, subtrial$immune.data)
    state$stageI.dose <- c(state$stageI.dose, subtrial$dose.data)
    return(state)
  }

  single_trial <- function(trial) {
    # cat("trial:", trial, "\n")
    state <- list(
      totaln = 0,
      npts = matrix(0, JJ, KK),
      ntox = matrix(0, JJ, KK),
      neff = matrix(0, JJ, KK),
      elimi = matrix(0, JJ, KK),
      total.duration = 0,
      stageI.pfs = NULL,
      stageI.immune = NULL,
      stageI.dose = NULL,
      npts.stage2 = matrix(0, JJ, KK),
      neff.stage2 = matrix(0, JJ, KK)
    )


    obdc <- cbind(selectdoseA = NA, selectdoseB = NA)

    startdose <- init.startdose[2]
    dosespace <- c(((1:(KK - 1)) - 1) * JJ + 1, (JJ * (KK - 1) + 1):(JJ * KK))
    subtriali <- 1

    ## ============ Stage I ============
    while (state$totaln < sum(ncohort) * cohortsize && subtriali <= JJ) {
      ret <- run_subtrial(dosespace, startdose, subtriali, state)
      subtrial <- ret$subtrial
      state <- update_state(state, subtrial)
      selectdose <- ret$selectdose # 子试验候选MTDC

      di <- one2two(selectdose, JJ, KK)[1]
      dj <- one2two(selectdose, JJ, KK)[2]

      # 第一次子试验的selectdose超过拐角处
      if ((subtriali == 1) &&
        (selectdose > JJ * (KK - 1) + 1) &&
        (selectdose < JJ * KK)) {
        startdose <- dj - 1
        dosespace <- -JJ + di + (1:(KK - 1) * JJ)
        subtriali <- subtriali + 1

        ret <- run_subtrial(dosespace, startdose, subtriali, state)
        subtrial <- ret$subtrial
        state <- update_state(state, subtrial)
        selectdose <- ret$selectdose
      }

      # 标记淘汰剂量（毒性过高）
      elimi.index <- dosespace[state$elimi[dosespace] == 1]
      for (x in elimi.index) {
        idx <- one2two(x, JJ, KK)
        state$elimi[idx[1], idx[2]] <- 1
      }


      if (selectdose == 99) {
        break
      } else {
        for (x in dosespace[dosespace > selectdose]) {
          idx <- one2two(x, JJ, KK)
          state$elimi[idx[1], idx[2]] <- 1
        }
      }
      # cat("state$elimi:\n")
      # print(state$elimi)
      # 下一个subtrial的起始剂量
      startdose <- ifelse(dj > 1, dj - 1, dj)
      dosespace <- -(JJ - 1) + subtriali + (1:(KK - 1) * JJ)
      subtriali <- subtriali + 1
    }
    npts <- state$npts
    ntox <- state$ntox
    neff <- state$neff
    elimi <- state$elimi
    total.duration <- state$total.duration
    stageI.dose <- state$stageI.dose
    stageI.pfs <- state$stageI.pfs
    stageI.immune <- state$stageI.immune
    npts.stage2 <- state$npts.stage2
    neff.stage2 <- state$neff.stage2

    # 剂量矩阵毒性概率估计
    phatT <- (ntox + 0.05) / (npts + 0.1) # 毒性概率后验均值
    phatT[elimi == 1 | npts == 0] <- 1.1
    phatT <- biviso(phatT, npts + 0.1, warn = TRUE)
    phatT <- matrix(as.numeric(phatT), JJ, KK)
    elimi[phatT > target] <- 1

    # 剂量矩阵疗效概率估计：P(pE < lower_e | data)>=cutoff.eli则淘汰
    for (j in 1:JJ) {
      for (k in 1:KK) {
        if (npts[j, k] > 0) {
          # cat("prob.tox.unsafe:",1-pbeta(target, ntox[j, k]+1, npts[j,k]-ntox[j,k]+1),"\n")
          prob.eff.unsafe <- pbeta(lower_e, neff[j, k] + 0.3, npts[j, k] - neff[j, k] + 0.7)
          # cat(j,k,"prob.eff.unsafe:", prob.eff.unsafe, "\n")
          if (prob.eff.unsafe >= cutoff.eli) {
            # cat(prob.eff.unsafe,"淘汰",j,k,"\n")
            elimi[j, k] <- 1
          }
        }
      }
    }
    adm.set <- (npts != 0) & (elimi == 0)
    adm.index <- which(adm.set == TRUE)

    # cat("npts:", "\n")
    # print(npts)
    # cat("ntox:", "\n")
    # print(ntox)
    # cat("neff:", "\n")
    # print(neff)
    # cat("elimi:", "\n")
    # print(elimi)
    #
    # cat("adm.index:", adm.index, "\n")
    ## ============ Stage II ============
    # cat("========Stage II Interim 1========\n")
    if (length(adm.index) == 0) {
      obdc.index <- 99
    } else {
      filter.index <- which(stageI.dose %in% adm.index)
      stageII.pfs <- stageI.pfs[filter.index]
      stageII.immune <- stageI.immune[filter.index]
      stageII.dose <- stageI.dose[filter.index]

      dose.table <- table(stageII.dose)

      for (i in seq_along(dose.table)) {
        di <- as.numeric(names(dose.table)[i])
        # Interim 1前每个剂量补齐患者至dose.max.sample/3
        addn <- dose.max.sample / 3 - dose.table[i]
        # cat("dose:", di, "addn:", addn, "\n")
        if (addn > 0) {
          stageII.pfs <- c(stageII.pfs, rlomax(addn, alpha.true, gamma.true[di]))
          stageII.yE.outcome <- rbinom(addn, 1, pE.true[di])
          stageII.immune <- c(stageII.immune, stageII.yE.outcome)
          stageII.dose <- c(stageII.dose, rep(di, addn))

          # 筛选剂量等于di的index并打印对应的pfs和immune
          # idx <- which(stageII.dose == di)
          # if (length(idx) > 0) {
          #   cat("median pfs for dose==", di,":", median(stageII.pfs[idx]), "\n")
          #   cat("sum immune for dose==", di,":", sum(stageII.immune[idx]), "\n")
          # }

          idx <- one2two(di, JJ, KK)
          di.x <- idx[1]
          di.y <- idx[2]
          npts.stage2[di.x, di.y] <- npts.stage2[di.x, di.y] + addn
          neff.stage2[di.x, di.y] <- neff.stage2[di.x, di.y] + sum(stageII.yE.outcome)
        }
      }

      # StageII Interim 1
      calc.res.interim1 <- calc.rmst.eff(
        pfs.data = stageII.pfs,
        immune.data = stageII.immune,
        dose.data = stageII.dose,
        JJ = JJ, KK = KK,
        tau = tau
      )
      rmst.interim1 <- calc.res.interim1$rmst
      pE.interim1 <- calc.res.interim1$pE
      # cat("stageII.pfs:",round(stageII.pfs,3),"\n")
      # cat("stageII.dose",stageII.dose,"\n")
      # cat("Stage II Interim 1 rmst:", "\n")
      # print(rmst.interim1)
      # cat("Stage II Interim 1 pE:", "\n")
      # print(pE.interim1)

      rmst.adm.interim1 <- rmst.interim1[adm.set]
      pE.adm.interim1 <- pE.interim1[adm.set]
      # 淘汰StageII Interim 1疗效过低和RMST过低的剂量组合
      adm.index.interim2 <- adm.index[which(
        pE.adm.interim1 >= lower_e & rmst.adm.interim1 >= rmst.min
      )]
      rmst.adm.interim2 <- rmst.adm.interim1[which(
        pE.adm.interim1 >= lower_e & rmst.adm.interim1 >= rmst.min
      )]
      # cat("adm.index.interim2:", "\n")
      # print(adm.index.interim2)

      # 限制可接受集合中RMST前K大的剂量进入Interim 2，K=(JJ*KK)/3
      if (length(adm.index.interim2) > 0 && length(adm.index.interim2) > floor((JJ * KK) / 3)) {
        K <- floor((JJ * KK) / 3)
        topK_idx <- order(rmst.adm.interim2, decreasing = TRUE)[1:K]
        adm.index.interim2 <- adm.index.interim2[topK_idx]
        rmst.adm.interim2 <- rmst.adm.interim2[topK_idx]
      }

      # StageII Interim 2
      # cat("========Stage II Interim 2========\n")
      if (length(adm.index.interim2) == 0) {
        obdc.index <- 99
      } else {
        dose.table2 <- table(factor(stageII.dose, levels = adm.index.interim2))
        for (i in seq_along(adm.index.interim2)) {
          di <- adm.index.interim2[i]
          # Interim 2前每个剂量补齐患者至2*dose.max.sample/3
          addn2 <- 2 * dose.max.sample / 3 - as.integer(dose.table2[i])
          # cat("dose:", di, "addn2:", addn2, "\n")
          if (addn2 <= 0) next

          stageII.pfs <- c(stageII.pfs, rlomax(addn2, alpha.true, gamma.true[di]))
          stageII.yE2 <- rbinom(addn2, 1, pE.true[di])
          stageII.immune <- c(stageII.immune, stageII.yE2)
          stageII.dose <- c(stageII.dose, rep(di, addn2))
          # cat("pfs:", stageII.pfs, "\n")
          # cat("immune:", stageII.immune, "\n")
          # cat("dose:", stageII.dose, "\n")
          
          # 筛选剂量等于di的index并打印对应的pfs和immune
          # idx <- which(stageII.dose == di)
          # if (length(idx) > 0) {
          #   cat("median pfs for dose==", di,":", median(stageII.pfs[idx]), "\n")
          #   cat("sum immune for dose==", di,":", sum(stageII.immune[idx]), "\n")
          # }

          idx <- one2two(di, JJ, KK)
          npts.stage2[idx[1], idx[2]] <- npts.stage2[idx[1], idx[2]] + addn2
          neff.stage2[idx[1], idx[2]] <- neff.stage2[idx[1], idx[2]] + sum(stageII.yE2)
        }

        # StageII Interim 2 分析
        calc.res.interim2 <- calc.rmst.eff(
          pfs.data = stageII.pfs,
          immune.data = stageII.immune,
          dose.data = stageII.dose,
          JJ = JJ, KK = KK,
          tau = tau
        )
        rmst.interim2 <- calc.res.interim2$rmst
        pE.interim2 <- calc.res.interim2$pE

        adm.set.interim2 <- matrix(FALSE, JJ, KK)
        adm.set.interim2[adm.index.interim2] <- TRUE
        rmst.adm.interim2 <- rmst.interim2[adm.set.interim2]
        pE.adm.interim2 <- pE.interim2[adm.set.interim2]

        # 淘汰Interim 2疗效过低的剂量
        adm.index.final <- adm.index.interim2[which(pE.adm.interim2 >= lower_e)]
        rmst.adm.final <- rmst.adm.interim2[which(pE.adm.interim2 >= lower_e)]

        # cat("stageII.pfs:",round(stageII.pfs,3),"\n")
        # cat("stageII.dose",stageII.dose,"\n")
        # cat("Stage II Interim 2 rmst:", "\n")
        # print(rmst.interim2)
        # cat("Stage II Interim 2 pE:", "\n")
        # print(pE.interim2)

        if (length(adm.index.final) == 0) {
          obdc.index <- 99
        } else {
          # Final Analysis: 补齐患者至dose.max.sample
          # cat("========Stage II Final========\n")
          dose.table3 <- table(factor(stageII.dose, levels = adm.index.final))
          for (i in seq_along(adm.index.final)) {
            di <- adm.index.final[i]
            addn3 <- dose.max.sample - as.integer(dose.table3[i])
            # cat("dose:", di, "addn3:", addn3, "\n")
            if (addn3 <= 0) next

            stageII.pfs <- c(stageII.pfs, rlomax(addn3, alpha.true, gamma.true[di]))
            stageII.yE3 <- rbinom(addn3, 1, pE.true[di])
            stageII.immune <- c(stageII.immune, stageII.yE3)
            stageII.dose <- c(stageII.dose, rep(di, addn3))

          # 筛选剂量等于di的index并打印对应的pfs和immune
          # idx <- which(stageII.dose == di)
          # if (length(idx) > 0) {
          #   cat("median pfs for dose==", di,":", median(stageII.pfs[idx]), "\n")
          #   cat("sum immune for dose==", di,":", sum(stageII.immune[idx]), "\n")
          # }

            idx <- one2two(di, JJ, KK)
            npts.stage2[idx[1], idx[2]] <- npts.stage2[idx[1], idx[2]] + addn3
            neff.stage2[idx[1], idx[2]] <- neff.stage2[idx[1], idx[2]] + sum(stageII.yE3)
          }

          # Final Stage II analysis
          calc.res <- calc.rmst.eff(
            pfs.data = stageII.pfs,
            immune.data = stageII.immune,
            dose.data = stageII.dose,
            JJ = JJ, KK = KK,
            tau = tau
          )
          rmst <- calc.res$rmst
          pE <- calc.res$pE

          adm.set.final <- matrix(FALSE, JJ, KK)
          adm.set.final[adm.index.final] <- TRUE
          rmst.adm.final <- rmst[adm.set.final]
          pE.adm.final <- pE[adm.set.final]

          # cat("stageII.pfs:",round(stageII.pfs,3),"\n")
          # cat("stageII.dose",stageII.dose,"\n")
          # cat("Stage II Final rmst:", "\n")
          # print(rmst)
          # cat("Stage II Final pE:", "\n")
          # print(pE)

          adm.index.update <- adm.index.final[which(pE.adm.final >= lower_e)]
          rmst.adm.update <- rmst.adm.final[which(pE.adm.final >= lower_e)]

          if (length(adm.index.update) == 0) {
            obdc.index <- 99
          } else {
            obdc.index <- adm.index.update[rmst.adm.update == max(rmst.adm.update)]
          }
        }
      }
      # cat("obdc.index:", obdc.index, "\n")
    }
    if (obdc.index == 99) {
      obdc <- matrix(c(99, 99), 1)
    } else {
      obdc <- matrix(one2two(obdc.index, JJ, KK), 1)
    }


    list(
      obdc = cbind(trial, obdc),
      nt = cbind(trial, npts),
      yt = cbind(trial, ntox),
      ye = cbind(trial, neff),
      duration = total.duration,
      nt2 = cbind(trial, npts.stage2),
      ye2 = cbind(trial, neff.stage2)
    )
  }

  plan(multisession, workers = max(1, detectCores() - 1)) # 并行
  # plan(multisession, workers = max(1, 20)) # 并行
  # plan(sequential) # 串行
  options(progressr.enable = TRUE)
  handlers(global = TRUE)
  handlers("txtprogressbar")

  res <- with_progress({
    p <- progressor(along = seq_len(ntrial))
    future_lapply(seq_len(ntrial), function(trial) {
      out <- single_trial(trial)
      p()
      return(out)
    },
    future.seed = TRUE,
    future.packages = "rjags"
    )
  })

  ntrial.obdc <- do.call(rbind, lapply(res, function(x) x[["obdc"]]))
  ntrial.nt <- do.call(rbind, lapply(res, function(x) x[["nt"]]))
  ntrial.yt <- do.call(rbind, lapply(res, function(x) x[["yt"]]))
  ntrial.ye <- do.call(rbind, lapply(res, function(x) x[["ye"]]))
  ntrial.duration <- sapply(res, function(x) x[["duration"]])
  ntrial.nt2 <- do.call(rbind, lapply(res, function(x) x[["nt2"]]))
  ntrial.ye2 <- do.call(rbind, lapply(res, function(x) x[["ye2"]]))

  ntrial.obdc <- data.frame(ntrial.obdc) # 存储着最终的OBDC
  obdc.selpercent <- matrix(0, nrow = JJ, ncol = KK)
  patpercent <- matrix(0, nrow = JJ, ncol = KK)

  # 将二维剂量组合转为一维编码，统计不同剂量组合被选为 OBDC 的频率
  for (triali in 1:ntrial) {
    obdcdata <- as.matrix(ntrial.obdc[ntrial.obdc$trial == triali, 2:3])
    if (obdcdata[1, 2] > KK) next
    obdclevel <- two2one(obdcdata[1, 1], obdcdata[1, 2], JJ)
    obdc.selpercent[obdclevel] <- obdc.selpercent[obdclevel] + 1
  }

  nptsdose <- matrix(0, nrow = JJ, ncol = KK)
  ntoxdose <- matrix(0, nrow = JJ, ncol = KK)
  neffdose <- matrix(0, nrow = JJ, ncol = KK)
  nptsdose2 <- matrix(0, nrow = JJ, ncol = KK)
  neffdose2 <- matrix(0, nrow = JJ, ncol = KK)

  for (i in seq(1, nrow(ntrial.yt), by = JJ)) {
    idx <- i + 0:(JJ - 1)
    # 提取第idx行，排除第1列
    nptsdose <- nptsdose + ntrial.nt[idx, -1]
    ntoxdose <- ntoxdose + ntrial.yt[idx, -1]
    neffdose <- neffdose + ntrial.ye[idx, -1]
    nptsdose2 <- nptsdose2 + ntrial.nt2[idx, -1]
    neffdose2 <- neffdose2 + ntrial.ye2[idx, -1]
  }

  out <- list(
    selpercent = format.matrix(obdc.selpercent / ntrial * 100, digits = 2),
    patpercent = format.matrix(nptsdose / sum(nptsdose) * 100, digits = 2),
    patpercent2 = format.matrix(nptsdose2 / sum(nptsdose2) * 100, digits = 2),
    patpercentall = format.matrix(
      (nptsdose + nptsdose2) / sum(nptsdose + nptsdose2) * 100,
      digits = 2
    ),
    npatients = format.matrix(nptsdose / ntrial, digits = 2),
    ntox = format.matrix(ntoxdose / ntrial, digits = 2),
    neff = format.matrix(neffdose / ntrial, digits = 2),
    npatients2 = format.matrix(nptsdose2 / ntrial, digits = 2),
    neff2 = format.matrix(neffdose2 / ntrial, digits = 2),
    totaln = format.num(sum(nptsdose / ntrial), digits = 1),
    totaltox = format.num(sum(ntoxdose / ntrial), digits = 1),
    totaleff = format.num(sum(neffdose / ntrial), digits = 1),
    totaln2 = format.num(sum(nptsdose2 / ntrial), digits = 1),
    totaleff2 = format.num(sum(neffdose2 / ntrial), digits = 1),
    duration = format.num(sum(ntrial.duration / ntrial), digits = 2)
  )
  return(out)
}
