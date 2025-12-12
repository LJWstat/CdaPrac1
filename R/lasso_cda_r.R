ST = function(z, lambda)
{
  if (z > lambda) return(z - lambda)
  if (z < -lambda) return(z + lambda)
  return(0)
}

lasso_cda_r = function(X, y, lambda, max_iter = 1000, tol = 1e-6) {
  n = nrow(X)
  p = ncol(X)
  beta = rep(0, p)
  r = y  # 잔차 벡터
  beta_old = beta

  for (iter in 1:max_iter)
  {
    max_change = 0

    for (j in 1:p)
    {
      X_j = X[, j]  # X의 j번째 열 추출
      X_j_norm2 = sum(X_j^2)  # ||X_j||^2 계산
      rho = sum(X_j * r) + X_j_norm2 * beta[j]

      beta_old[j] = beta[j]
      beta[j] = ST(rho / X_j_norm2, lambda / X_j_norm2)

      r = r + X_j * (beta_old[j] - beta[j])  # 잔차 업데이트
      max_change = max(max_change, abs(beta[j] - beta_old[j]))
    }

    if (max_change < tol) break
  }

  return(beta)
}
