#' Soft-thresholding Operator
#'
#' Soft-thresholding 함수는 Lasso의 좌표강하(Coordinate Descent Algorithm)에서
#' 회귀계수 업데이트에 사용되는 핵심 연산입니다.
#'
#' @param z Numeric scalar. 업데이트 전의 좌표 값.
#' @param lambda Numeric scalar. L1 penalty 계수.
#'
#' @return Numeric scalar. Soft-thresholding 결과값.
#'
#' @export
#'
#' @examples
#' ST(3, 1)
#' ST(-0.5, 1)
#' ST(0.2, 0.5)
ST = function(z, lambda)
{
  if (z > lambda) return(z - lambda)
  if (z < -lambda) return(z + lambda)
  return(0)
}

#' Lasso Regression via Coordinate Descent
#'
#' 이 함수는 L1-penalized 회귀(Lasso)를 좌표강하법(Coordinate Descent Algorithm, CDA)
#' 으로 구현합니다. 각 좌표(변수)를 순차적으로 업데이트하며, soft-thresholding을
#' 이용해 회귀계수를 효율적으로 추정합니다.
#'
#' @param X Numeric matrix (n × p). 디자인 행렬.
#' @param y Numeric vector (n). 반응변수.
#' @param lambda Numeric scalar. L1 penalty 계수.
#' @param max_iter Maximum number of coordinate descent iterations.
#' @param tol Convergence tolerance. 계수 변화량이 tol 미만이면 수렴으로 판단.
#'
#' @return Numeric vector of estimated coefficients (length p).
#'
#' @details
#' Lasso 추정량은 다음 최적화 문제를 푸는 것과 동일합니다:
#'
#' \deqn{\hat{\beta} = \arg\min_\beta \frac{1}{2}\|y - X\beta\|_2^2 + \lambda\|\beta\|_1}
#'
#' 좌표강하법에서는 한 번에 하나의 변수 \eqn{\beta_j} 만 업데이트하며,
#' soft-thresholding 연산자가 사용됩니다.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(100 * 10), 100, 10)
#' beta_true <- c(3, -2, rep(0, 8))
#' y <- X %*% beta_true + rnorm(100)
#'
#' lasso_cda_r(X, y, lambda = 0.5)
#'
#' @export
lasso_cda_r = function(X, y, lambda, max_iter = 1000, tol = 1e-6) {
  n = nrow(X)                     # 관측치 개수 n (행 개수) 계산
  p = ncol(X)                     # 변수 개수 p (열 개수) 계산
  beta = rep(0, p)                # 회귀계수 초기값을 0 벡터로 설정 (길이 p)
  r = y                           # 초기 잔차 r = y - X %*% beta 이지만 beta=0이므로 r = y
  beta_old = beta                 # 이전 단계의 계수를 저장할 벡터 초기화

  for (iter in 1:max_iter)        # 최대 max_iter 번까지 좌표강하 반복
  {
    max_change = 0                # 이번 iteration에서 계수 변화량의 최대값을 추적

    for (j in 1:p)                # 각 좌표(j번째 변수)에 대해 순차적으로 업데이트
    {
      X_j = X[, j]                # 디자인 행렬 X의 j번째 열 (변수 j) 추출
      X_j_norm2 = sum(X_j^2)      # ||X_j||^2 = sum(X_j^2), 좌표강하 공식에 사용됨
      rho = sum(X_j * r) + X_j_norm2 * beta[j]
      # rho = X_j^T (잔차 + X_j * beta_j) = X_j^T (y - X_{-j} beta_{-j})

      beta_old[j] = beta[j]       # 업데이트 전의 beta_j 값을 저장
      beta[j] = ST(rho / X_j_norm2, lambda / X_j_norm2)
      # soft-thresholding 적용해 새로운 beta_j 계산

      r = r + X_j * (beta_old[j] - beta[j])
      # 잔차 r = y - X beta 를 효율적으로 업데이트
      # (직접 y - X%*%beta 다시 계산하는 대신, j좌표 변화분만 반영)
      max_change = max(max_change, abs(beta[j] - beta_old[j]))
      # 이번 좌표 업데이트에서의 변화량 |beta_j - beta_j_old| 중 최대값 갱신
    }

    if (max_change < tol) break   # 모든 좌표 변화량이 tol보다 작으면 수렴했다고 보고 종료
  }

  return(beta)                    # 최종 추정된 Lasso 회귀계수 벡터 반환
}
