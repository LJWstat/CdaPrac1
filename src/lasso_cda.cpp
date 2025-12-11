#include <Rcpp.h>                       // Rcpp 기능 사용
using namespace Rcpp;                   // Rcpp 네임스페이스 사용

// -------------------------------
// Soft-thresholding 함수
// -------------------------------

// [[Rcpp::export]]
double soft_threshold(double z, double lambda) {
  if (z > lambda) return z - lambda;    // z가 λ보다 크면 z - λ
  if (z < -lambda) return z + lambda;   // z가 -λ보다 작으면 z + λ
  return 0.0;                           // |z| ≤ λ 이면 0 (Lasso의 핵심 shrinkage)
}


// -------------------------------
// Coordinate Descent Lasso 알고리즘
// -------------------------------

// [[Rcpp::export]]
NumericVector lasso_cda_cpp(NumericMatrix X, NumericVector y,
                            double lambda, int max_iter = 1000, double tol = 1e-6)
{
  int n = X.nrow(),                    // 행 개수: 샘플 수 n
    p = X.ncol();                    // 열 개수: 변수 수 p

  NumericVector beta(p, 0.0);          // 회귀계수 β 초기값 (모두 0)
  NumericVector X_j(n);                // X의 j번째 열을 저장할 임시 벡터
  NumericVector r = clone(y);          // 초기 잔차 r = y   (β = 0이므로)

  double beta_old;                     // 업데이트 전 β_j 저장용
  double max_change;                   // 한 iteration에서의 최대 계수 변화량

  // -------------------------------
  // 좌표강하 반복 시작
  // -------------------------------
  for (int iter = 0; iter < max_iter; iter++) {
    max_change = 0.0;

    // ----------------------------
    // 각 좌표 j = 1..p 업데이트
    // ----------------------------
    for (int j = 0; j < p; j++)
    {
      // X_j ← X[, j] (j번째 변수)
      for (int i = 0; i < n; i++) {
        X_j[i] = X(i, j);              // j열을 복사해 column vector X_j 구성
      }

      // ||X_j||^2 = sum(X_j^2)
      double X_j_norm2 = sum(X_j * X_j);

      // ρ_j = X_j^T r + ||X_j||^2 * β_j
      // (partial residual 기반의 closed-form 업데이트)
      double rho = sum(X_j * r) + X_j_norm2 * beta[j];

      beta_old = beta[j];              // 업데이트 이전 β_j 저장

      // soft-threshold 적용하여 β_j 업데이트
      beta[j] = soft_threshold(rho / X_j_norm2,
                               lambda / X_j_norm2);

      // 잔차 업데이트: r ← r + X_j * (β_old - β_new)
      // 전체 r = y - Xβ 를 다시 계산하지 않고 효율적으로 업데이트
      r = r + X_j * (beta_old - beta[j]);

      // 이번 좌표 업데이트의 변화량 기록
      max_change = std::max(max_change,
                            std::abs(beta[j] - beta_old));
    }

    // ----------------------------
    // 수렴 검사: 한 번의 전체 순회에서 변화 ≤ tol
    // ----------------------------
    if (max_change < tol) break;
  }

  return beta;                         // 최종 Lasso 계수 반환
}
