#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double soft_threshold(double z, double lambda) {
  if (z > lambda) return z - lambda;
  if (z < -lambda) return z + lambda;
  return 0.0;
}

// [[Rcpp::export]]
NumericVector lasso_cda_cpp(NumericMatrix X, NumericVector y, double lambda, int max_iter = 1000, double tol = 1e-6) {
  int n = X.nrow(), p = X.ncol();
  NumericVector beta(p, 0.0);
  NumericVector X_j(n);
  NumericVector r = clone(y);
  double beta_old, max_change;

  for (int iter = 0; iter < max_iter; iter++) {
    max_change = 0.0;

    for (int j = 0; j < p; j++)
    {
      for (int i = 0; i < n; i++)
      {
        X_j[i] = X(i, j);
      }

      double X_j_norm2 = sum(X_j * X_j);
      double rho = sum(X_j * r) + X_j_norm2 * beta[j];
      beta_old = beta[j];

      beta[j] = soft_threshold(rho / X_j_norm2, lambda / X_j_norm2);

      r = r + X_j * (beta_old - beta[j]);
      max_change = std::max(max_change, std::abs(beta[j] - beta_old));
    }

    if (max_change < tol) break;
  }
  return beta;
}
