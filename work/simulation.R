set.seed(42)

# 1. 시뮬레이션 설정 ----
n = 100   # 샘플 개수
p = 10    # 변수 개수
sparsity = 5  # 실제로 신호를 갖는 변수 개수
sigma = 1  # 노이즈 크기
lambda = 8  # Lasso 페널티 파라미터

# 2. 데이터 생성 ----
X = matrix(rnorm(n * p), n, p)  # X ~ N(0,1)
beta_true = rep(0, p)            # 초기값 (모든 계수를 0으로 설정)
beta_true[sample(1:p, sparsity)] = runif(sparsity, 1, 3)  # 일부 변수만 신호 부여
beta_0 = 2
y = beta_0 + X %*% beta_true + rnorm(n, 0, sigma)  # y 생성

# 3. Lasso 실행 ----
lasso_result = lasso_cda(X, y, lambda)

# 4. 결과 평가 ----
beta_hat = lasso_result$beta  # 추정된 회귀 계수
intercept_hat = lasso_result$intercept  # 추정된 절편

# Mean Squared Error (MSE) 계산
y_pred = X %*% beta_hat + intercept_hat
mse = mean((y - y_pred)^2)

# R-squared 계산
r_squared = 1 - sum((y - y_pred)^2) / sum((y - mean(y))^2)

# 선택된 변수 평가 (True Positive, False Positive)
selected_vars = which(beta_hat != 0)
true_vars = which(beta_true != 0)

TP = length(intersect(selected_vars, true_vars))  # 실제 신호 변수 중 선택된 변수
FP = length(setdiff(selected_vars, true_vars))    # 신호가 없는데 선택된 변수
FN = length(setdiff(true_vars, selected_vars))    # 신호가 있는데 선택되지 않은 변수

# 5. 결과 출력 ----
cat("\n True Beta:\n", beta_true, "\n")
cat("\n Estimated Beta:\n", beta_hat, "\n")
cat("\n Intercept:\n", intercept_hat, "\n")
cat("\n MSE:", mse, "\n")
cat("R-squared:", r_squared, "\n")
cat("TP:", TP, "| FP:", FP, "| FN:", FN, "\n")
