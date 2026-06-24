setwd("R")
source("get.oc.comb.R")
source("calc.lomax.scale.R")
source("scenarios.R")
# source("scenariosAB.R")
seed <- 10
set.seed(seed)

scenario <- 1
tite <- FALSE
# tite <- TRUE

current_time <- format(Sys.time(), "%Y%m%d_%H%M%S")
# current_time <- "20260622_sample24"

for (start_high in c(FALSE, TRUE)) {
  # start_high <- FALSE
  for (scenario in c(1:12)) {
    for (dose.max.sample in c(36)) {
    pT.true <- list.p.true.tox[[scenario]]
    pE.true <- list.p.true.eff[[scenario]]
    rmst.true <- list.rmst.true[[scenario]]
    true.obdc <- list.tdose[[scenario]]

    startdose <- if (start_high) {
      c(1, ncol(list.p.true.tox[[scenario]]))
    } else {
      c(1, 1)
    }

    cat("Scenario:", scenario, "\n")
    cat("Startdose:", paste0("(", startdose[1], ",", startdose[2], ")"), "\n")
    cat("TITE:", tite, "\n")
    cat("dose.max.sample:", dose.max.sample, "\n")
    alpha.true <- 1.5
    gamma.true <- matrix(NA, nrow = nrow(rmst.true), ncol = ncol(rmst.true))
    tau <- 12
    for (i in 1:nrow(rmst.true)) {
      for (j in 1:ncol(rmst.true)) {
        gamma.true[i, j] <- calc.lomax.scale(alpha.true, tau, rmst.true[i, j])
      }
    }

    target <- 0.3
    lower_e <- 0.4
    rmst.min <- 4
    dose.ncohort <- 2
    # 每个子试验cohort数量
    ncohort.first.subtrial <- dose.ncohort * (ncol(pT.true) + nrow(pT.true) - 1)
    ncohort.other.subtrial <- dose.ncohort * (ncol(pT.true) - 1)
    ncohort <- c(
      ncohort.first.subtrial,
      rep(ncohort.other.subtrial, nrow(pT.true) - 1)
    )
    cohortsize <- 3
    # n.earlystop <- 6
    n.earlystop <- 9

    ntrial <- 1000


    start.flag <- ifelse(startdose[2] == 1, ".low", ".high")
    output.dir <- paste0("../output/", current_time, "_d", dose.max.sample, "_es", n.earlystop, "/")
    if (!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)

    start_time <- Sys.time()

    oc.comb <- get.oc.comb(
      target = target, lower_e = lower_e, rmst.min = rmst.min, pT.true = pT.true,
      pE.true = pE.true, rmst.true = rmst.true, alpha.true = alpha.true,
      gamma.true = gamma.true, ncohort = ncohort, cohortsize = cohortsize,
      n.earlystop = n.earlystop, startdose = startdose, tau = tau,
      dose.max.sample = dose.max.sample, ntrial = ntrial, tite = tite
    )

    filename <- scenario
    if (tite) filename <- paste0(filename, ".tite")
    filename <- paste0(filename, start.flag)
    full.filepath <- paste(output.dir, filename, ".txt", sep = "")
    cat("Output file:", full.filepath, "\n")
    sink(full.filepath, append = F)
    cat("startdose:", paste0("(", startdose[1], ",", startdose[2], ")"), "\n")
    cat("target:", target, "\n")
    cat("ncohort:", ncohort, "\n")
    cat("cohortsize:", cohortsize, "\n")
    cat("Stage I max smaple size:", sum(ncohort) * cohortsize, "\n")
    cat("Stage II max sample size for each dose:", dose.max.sample, "\n")
    cat("n.earlystop:", n.earlystop, "\n")
    cat("ntrial:", ntrial, "\n")
    cat("tite:", tite, "\n")
    cat("\n")

    cat("scenario:", scenario, "\n")
    cat("true toxicity probability:\n")
    print(pT.true)
    cat("true efficacy probability:\n")
    print(pE.true)
    cat("true RMST:\n")
    print(rmst.true)
    cat("true obdc:", true.obdc, "\n")
    cat("\n")

    cat("selection percentage at each dose combination (%):\n")
    print(oc.comb$selpercent)
    cat("select OBDC correctly (%):\n")
    if (true.obdc[1] != -1) {
      print(sum(oc.comb$selpercent[true.obdc]))
    } else {
      print("null")
    }
    cat("early stop percentage (%):\n")
    print(100 - sum(oc.comb$selpercent))
    cat("\n")

    cat("average percent of patients(%):\n")
    print(oc.comb$patpercentall)
    cat("average percent of patients treated at OBDC correctly(%):\n")
    if (true.obdc[1] != -1) {
      print(sum(oc.comb$patpercentall[true.obdc]))
    } else {
      print("null")
    }


    cat("Stage I average number of patients:\n")
    print(oc.comb$npatients)
    cat("Stage I average number of patients treated at OBDC correctly:\n")
    if (true.obdc[1] != -1) {
      print(sum(oc.comb$npatients[true.obdc]))
    } else {
      print("null")
    }
    cat("Stage I total number of patients:\n")
    print(oc.comb$totaln)
    cat("\n")


    cat("Stage II average number of patients:\n")
    print(oc.comb$npatients2)
    cat("Stage II average number of patients treated at OBDC correctly:\n")
    if (true.obdc[1] != -1) {
      print(sum(oc.comb$npatients2[true.obdc]))
    } else {
      print("null")
    }
    cat("Stage II total number of patients:\n")
    print(oc.comb$totaln2)
    cat("\n")


    cat("Stage I average percent of patients(%):\n")
    print(oc.comb$patpercent)
    cat("Stage II average percent of patients(%):\n")
    print(oc.comb$patpercent2)
    cat("\n")

    cat("Stage I average number of toxicity:\n")
    print(oc.comb$ntox)
    cat("Stage I total number of toxicity:\n")
    print(oc.comb$totaltox)
    cat("\n")

    cat("Stage I average number of efficacy:\n")
    print(oc.comb$neff)
    cat("Stage I total number of efficacy:\n")
    print(oc.comb$totaleff)
    cat("\n")

    cat("Stage II average number of efficacy:\n")
    print(oc.comb$neff2)
    cat("Stage II total number of efficacy:\n")
    print(oc.comb$totaleff2)
    cat("\n")


    cat("average trial duration:\n")
    print(sum(oc.comb$duration))

    end_time <- Sys.time()
    print(end_time - start_time)
    sink()
    }
  }
}
