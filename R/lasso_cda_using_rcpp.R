
#'
#' @export
#'
lasso_cda = function(X, y, lambda, max_iter = 1000, tol = 1e-6)
{
  # 표준화
  X_scaled = scale(X, center = TRUE, scale = TRUE)
  y_mean = mean(y)
  y_scaled = y - y_mean

  # Rcpp를 이용해 coordinate descent 실행
  beta_standardized = lasso_cda_cpp(X_scaled, y_scaled, lambda, max_iter, tol)

  # 원래 스케일로 변환
  X_means = attr(X_scaled, "scaled:center")
  X_sds = attr(X_scaled, "scaled:scale")

  beta_original = beta_standardized / X_sds
  intercept = y_mean - sum(beta_original * X_means)

  return(list(beta = beta_original, intercept = intercept))
}
