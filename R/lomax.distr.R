# Lomax分布的概率密度函数
dlomax <- function(t, alpha, gamma) {
  d.val <- alpha * gamma^alpha / (t + gamma)^(alpha + 1)
  return(d.val)
}

# 使用逆变换抽样法从Lomax分布中抽取样本
rlomax <- function(n, alpha, gamma) {
  u <- runif(n, min = 0, max = 1)
  r.val <- gamma * ((1 - u)^(-1 / alpha) - 1)
  return(r.val)
}

# Lomax分布的生存函数
slomax <- function(t, alpha, gamma) {
  s.val <- (gamma / (t + gamma))^alpha
  return(s.val)
}
