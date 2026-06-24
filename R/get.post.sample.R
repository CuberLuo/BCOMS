get.post.sample <- function(pfs.data, immune.data, dose.data, n.pat, ndose, tau, silent = TRUE) {
  time <- pmin(pfs.data, tau) # 观测时间
  status <- as.numeric(pfs.data <= tau) # 1=事件发生，0=右删失

  dat.list <- list(
    immune.data = immune.data,
    dose.data = dose.data,
    n.pat = n.pat,
    ndose = ndose,
    zeros = rep(0, n.pat), # 伪观测，全是 0
    time = time,
    status = status
  )

  jags.inits <- function(chain) {
    list(
      alpha = rep(1, ndose),
      gamma = rep(4, ndose),
      beta = log(0.1) / 2,
      pE = rep(0.3, ndose),
      .RNG.name = ifelse(chain == 1, "base::Super-Duper",
        "base::Wichmann-Hill"
      ), # 随机数算法
      .RNG.seed = chain
    )
  }

  suppressWarnings({
    capture.output({ # 屏蔽输出Finished running the simulation
      mcmc.model <- run.jags(
        data = dat.list,
        inits = jags.inits,
        model = "survival.model.bugs",
        monitor = c("alpha", "gamma", "pE"), # "beta"暂时不需要后验
        n.chains = 2,
        sample = 2500,
        burnin = 1000,
        method = "rjags",
        silent.jags = silent
      )
    })
  })


  param.mcmc <- as.mcmc.list(mcmc.model)
  param.df <- do.call(rbind, param.mcmc)
  return(param.df)
}
