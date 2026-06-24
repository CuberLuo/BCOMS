calc.lomax.scale <- function(shape.val, tau, rmst) {
  alpha <- shape.val
  f <- function(gamma, tau, rmst) {
    if (alpha == 1) {
      res <- gamma * log(tau / gamma + 1)
    } else {
      res <- gamma / (1 - alpha) * ((1 + tau / gamma)^(1 - alpha) - 1) - rmst
    }
    return(res)
  }
  gamma <- nleqslv(x = 1, fn = f, tau = tau, rmst = rmst)$x
  return(gamma)
}
