gen.tite <- function(n, pi, dist = 1, alpha = 0.5, Tobs = 1) {
  weib <- function(n, pi, pihalft) {
    ## solve parameters for Weibull given pi=1-S(T) and phalft=1-S(T/2)
    alpha <- log(log(1 - pi) / log(1 - pihalft)) / log(2)
    lambda <- -log(1 - pi) / (Tobs^alpha)
    t <- (-log(runif(n)) / lambda)^(1 / alpha)
    return(t)
  }

  llogit <- function(n, pi, pihalft) {
    ## solve parameters for log-logistic given pi=1-S(T) and phalft=1-S(T/2)
    alpha <- log((1 / (1 - pi) - 1) / (1 / (1 - pihalft) - 1)) / log(2)
    lambda <- (1 / (1 - pi) - 1) / (Tobs^alpha)
    t <- ((1 / runif(n) - 1) / lambda)^(1 / alpha)
    return(t)
  }


  event <- rep(0, n)
  t.event <- rep(0, n)

  #### uniform
  if (dist == 1) { # 50% event in (0, 1/2T)
    event <- rbinom(n, 1, pi)
    nevent.st <- sum(event)
    t.event[event == 0] <- Tobs
    t.event[event == 1] <- runif(nevent.st, 0, Tobs)
  }
  #### Weibull
  if (dist == 2) {
    pihalft <- alpha * pi # alpha*100% event in (0, 1/2T)
    t.event <- weib(n, pi, pihalft)
    event[t.event <= Tobs] <- 1
    nevent.st <- sum(event)
    t.event[event == 0] <- Tobs
  }
  #### log-logistic
  if (dist == 3) {
    pihalft <- alpha * pi # alpha*100% event in (0, 1/2T)
    t.event <- llogit(n, pi, pihalft)
    event[t.event <= Tobs] <- 1
    nevent.st <- sum(event)
    t.event[event == 0] <- Tobs
  }
  return(list(event = event, t.event = t.event, nevent.st = nevent.st))
}
